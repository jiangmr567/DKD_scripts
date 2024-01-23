#########################################################################
# File Name: callpeak.R
# 添加peak矩阵，并富集motif
#########################################################################
# load library
library(ArchR)
library(parallel)
library(Seurat)
library(pryr) 
library(cowplot)
addArchRThreads(threads = 16)
set.seed(123)
library("BSgenome.Mfascicularis.NCBI.5.0")
assignInNamespace(".h5read", function (file = NULL, name = NULL, method = "fast", index = NULL, 
    start = NULL, block = NULL, count = NULL) {
    res <- rhdf5::h5read(file = file, name = name, index = index, 
            start = start, block = block, count = count)
    o <- h5closeAll()
    return(res)
}, ns="ArchR")
seqnames(BSgenome.Mfascicularis.NCBI.5.0)=gsub("MFA","chr",seqnames(BSgenome.Mfascicularis.NCBI.5.0))

args = commandArgs(T)
outdir<-args[1]
setwd(outdir)

genomeAnnotation <- readRDS("../data/Macaca_genomeAnnotation")
geneAnnotation <- readRDS("../data/Macaca_geneAnnotation")

proj <- readRDS("Save-ArchR-Project.rds")

system("pip install macs2")
pathToMacs2 <- findMacs2()

# peak calling with macs2
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "celltype",force = TRUE)

proj <- addReproduciblePeakSet(
    ArchRProj = proj, 
    groupBy = "celltype", 
    pathToMacs2 = pathToMacs2,
    cutOff=0.1,
    genomeSize=3e9,
    maxPeaks = 500000,
    force = TRUE
)
proj <- addPeakMatrix(proj)
allpeaks <- getPeakSet(proj)
df <- DataFrame(allpeaks)
head(df)

proj <- addImputeWeights(proj,seed=1)

# celltype specific genes (for manuscript)
markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "celltype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
dir.create("markerGene_celltype_final")
saveRDS(markersGS, file="markerGene_celltype_final/markerGene_celltype_final.rds")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")
printMKG <- function(MKGdata){
  outdata <- matrix(rep(0,nrow(MKGdata)*ncol(MKGdata)),nrow=nrow(MKGdata))
  colnames(outdata) <- colnames(MKGdata)
  for(i in 1:ncol(outdata)){
    outdata[,i] <- as.vector(MKGdata[,i])
  }
  return(outdata)
}

for(Group in names(markerList)){
  write.table(printMKG(markerList[Group][[1]]), file=paste0("markerGene_celltype_final/",Group,"_markerGene.txt"),row.names=F,col.names=T,sep="\t",quote=F)
}


# celltype specific peaks (for manuscript)
dir.create("markerPeak_celltype_final")
markersPeaks <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "PeakMatrix", 
    groupBy = "celltype",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, file="markerPeak_celltype_final/markerPeak_celltype_final.rds")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

output_diffpeak <- function(this_data,cellname){
    out_data <- cbind(as.vector(this_data@seqnames),
                      this_data@ranges@start,
                      this_data@ranges@start + this_data@ranges@width,
                      paste0(cellname,"_",seq(length(this_data@seqnames))),
                      this_data$Log2FC,
                      this_data$FDR)
    write.table(out_data,file=paste0("markerPeak_celltype_final/",cellname,"_markerPeak.bed"),quote=F,row.names=F,col.names=F,sep="\t")
}
for(i in names(markerList)){
	output_diffpeak(markerList[[i]],i)}
    
# motif annotation
motifs <- TFBSTools::getMatrixSet(x = JASPAR2020::JASPAR2020,
          opts = list(species = 9606, collection = "CORE"))
obj <- ArchR:::.summarizeJASPARMotifs(motifs)
proj <- addMotifAnnotations(
  ArchRProj = proj, 
  motifSet = "Custom", 
  name = paste(f,"Motif",sep='-'),
  motifPWMs = TFBSTools::toPWM(obj$motifs),
  force = TRUE)

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
saveRDS(enrichMotifs, file="markerPeak_celltype_final/celltype_enrichMotifs.rds")

saveRDS(proj, "Save-ArchR-Project.rds")

levels <- rev(c("Immune", "PODO", "ICB", "ICA", "Myofib", "Endo", "DCT", 
         "CNT", "PC", "PEC", "TAL","PT_VCAM1","PT"))
         
markersGS <- readRDS("./markerGene_celltype_final/markerGene_celltype_final.rds")
p1heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS[,levels], 
  cutOff = "FDR <= 0.01 & Log2FC >= 1", 
  #labelMarkers = markerGenes,
  nLabel=0,
  clusterCols = FALSE,
  transpose = TRUE
)
p1 <- ComplexHeatmap::draw(p1heatmapGS,heatmap_legend_side = "bot", annotation_legend_side = "bot")
p1
plotPDF(p1, name = "markerGene_celltype_GSheatmap_final.pdf", ArchRProj = proj, addDOC = FALSE, width = 4, height = 6)

# celltype specific peaks (for manuscript)
markersPeaks <- readRDS("./markerPeak_celltype_final/markerPeak_celltype_final.rds")
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks[,levels], 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  nLabel=0,clusterCols = FALSE,
  transpose = TRUE
)
p2 <- draw(heatmapPeaks,heatmap_legend_side = "bot", annotation_legend_side = "bot")
p2
plotPDF(p2, name = "markerPeak_celltype_SIGheatmap_final.pdf", width = 4, height = 6, ArchRProj = proj, addDOC = FALSE)

# Rename motifs for more aesthetic plotting:
rownames(enrichMotifs) <- lapply(rownames(enrichMotifs), function(x) strsplit(x, "_")[[1]][1]) %>% unlist()
# Subset to clusters that have at least some enrichment
log10pCut <- 10
levels <- rev(c("Immune", "PODO", "ICB", "ICA", "Myofib", "Endo", "DCT", 
                "CNT", "PC", "PEC", "TAL","PT_VCAM1","PT"))
#ArchR Heatmap
heatmapEM <- plotEnrichHeatmap(
  enrichMotifs[,levels], 
  n=5, clusterCols=FALSE,
  transpose=TRUE, 
  cutOff=10
)
pdf("enrichMotifs_celltype_Heatmap_final.pdf", width = 8, height = 5)
ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()