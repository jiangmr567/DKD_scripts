#########################################################################
# File Name: motif_chromVAR.R
# 计算所有motif注释中的每个细胞的偏差,以及关键TF的足迹分析
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

# chromVAR偏差富集
proj <- addBgdPeaks(proj)

proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "chromVAR_VariableMotifDeviationScores_final.pdf", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

# plot key motifs/TFs as example
motifs <- c("HNF4G","HNF4A","HNF1A","HNF1B","PPARD","RXRG","EHF","ETV1","ETV4","ELF1","IKZF1","GABPA")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% paste0("z:",c("RXRG.var.2_441","RARA..RXRG_349"))]
markerMotifs

proj <- addImputeWeights(proj,seed=1)
##MotifMatrix
p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAPHarmony",
    imputeWeights = getImputeWeights(proj)
)
plotPDF(p, name = "motif_chromVAR_UMAP.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

proj <- addGroupCoverages(ArchRProj = proj,minCells = 50,groupBy = "group",force = TRUE)
saveRDS(proj,"Save-ArchR-Project.rds")

###PT足迹分析
motifPositions <- getPositions(proj)
motifs <- c("HNF4G","HNF4A","HNF1A","HNF1B","PPARD","RXRG")
motifPositions <- getPositions(proj)
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% c("HNF4A.var.2_383","RARA..RXRG_349","RXRG.var.2_441")]
markerMotifs

peak_gr <- readRDS("./PeakCalls/PT-reproduciblePeaks.gr.rds")
peak_gr

for(i in markerMotifs){
    x=findOverlaps(motifPositions[[i]],peak_gr)
    d=as.data.frame(motifPositions[[i]][unique(queryHits(x))])
    print(i)
    print(dim(d))
    motifPositions[[i]] <- GRanges(d)
}

seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "group",
  #minCells = 50
)
seFoot<- seFoot[,grepl("PT-", colnames(seFoot))]
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "PT-Footprints-Subtract",
  addDOC = FALSE,
  smoothWindow = 5
)

##Endo足迹分析
motifPositions <- getPositions(proj)
motifs <- c("EHF","ETV1","ETV4","ELF1","IKZF1","GABPA")
motifPositions <- getPositions(proj)
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% c("HNF4A.var.2_383","THRB.var.2_460","THRB.var.3_461")]
markerMotifs
peak_gr <- readRDS("./PeakCalls/Endo-reproduciblePeaks.gr.rds")
for(i in markerMotifs){
    x=findOverlaps(motifPositions[[i]],peak_gr)
    d=as.data.frame(motifPositions[[i]][unique(queryHits(x))])
    print(i)
    print(dim(d))
    motifPositions[[i]] <- GRanges(d)
}
seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "group",
  #minCells = 50
)
seFoot<- seFoot[,grepl("Endo-", colnames(seFoot))]
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "Endo-Footprints-Subtract",
  addDOC = FALSE,
  smoothWindow = 5
)

##提取motif矩阵
motifMatrix <- getMatrixFromProject(proj, useMatrix="MotifMatrix")
markerMotifs <- c('HNF4G_583','HNF4A_582','HNF1B_86','ETV1_568','IKZF1_397','GABPA_576')
motif_df <- as.data.frame(t(assays(motifMatrix)$deviations[markerMotifs,]))
meta_df <- proj@cellColData[,c('state','celltype','group')]
head(meta_df)
motif_df$cellname <- rownames(motif_df)
meta_df$cellname <- rownames(meta_df)
plot_df <- merge(motif_df,meta_df,by='cellname')
plot_df <- as.data.frame(plot_df)[,c(-1,-8,-10)]
head(plot_df)
#使用aggregate函数取平均
result=aggregate(.~celltype,plot_df,mean)
rownames(result) <- result$celltype
result <- t(result[,-1])
result
#对数据按照基因scale
Result_matrix = t(scale(t(result)))
#Result_matrix = MinMax(Result_matrix, min = -4, max = 4)
levels=c("Immune", "PODO", "ICB", "ICA", "Myofib", "Endo", "DCT", 
         "CNT", "PC", "PEC", "TAL","PT_VCAM1","PT")
Result_matrix <- Result_matrix[,rev(levels)]
rownames(Result_matrix) <- gsub('[_].*','',rownames(Result_matrix))
Result_matrix

library(ComplexHeatmap)
p1 <- Heatmap(Result_matrix,name="Scores",
         col=c("#4264a3","white","#9b1a28"),
         cluster_rows = F,
         cluster_columns = F,
         show_column_names = T,
         show_row_names = T,
         column_title = "DeviationScores of Motif",
         row_names_side = "left",column_names_side = "bottom",
         column_names_max_height = unit(200, "cm"))
pdf("./Plots/motif_chromVAR_heatmap.pdf",width = 6,height = 3.5)
p1
dev.off()