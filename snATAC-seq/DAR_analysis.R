#########################################################################
# File Name: DAR_analysis.R
# 疾病差异peak分析
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

# disease v.s. normal DAR
dir.create("DAR_celltype_binomial")
diseaseMKG <- function(cellname){
  proj$tmpCluster <- "Other"
  proj@cellColData[proj$state=='DM' & proj$celltype== cellname,]$tmpCluster <- "DM"
  proj@cellColData[proj$state=='CTRL' & proj$celltype== cellname,]$tmpCluster <- "CTRL"
  print(cellname)
  print(table(proj$tmpCluster))

  ############# marker peak
  DMPeaks_tmp <- getMarkerFeatures(
      ArchRProj = proj, 
      useMatrix = "PeakMatrix", 
      groupBy = "tmpCluster",
    useGroups="DM",
    bgdGroups="CTRL",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod="binomial",binarize=TRUE
  )
  DMList_peak_tmp <- getMarkers(DMPeaks_tmp, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", returnGR = TRUE)
  this_data <- DMList_peak_tmp$DM
  out_data <- cbind(as.vector(this_data@seqnames),
                    this_data@ranges@start,
                    this_data@ranges@start + this_data@ranges@width,
                    paste0(cellname,"_",seq(length(this_data@seqnames))),
                    this_data$Log2FC,
                    this_data$FDR)
  write.table(out_data,file=paste0("DAR_celltype_binomial/",celltype,"_DAR.bed"),quote=F,row.names=F,col.names=F,sep="\t")
}
for(celltype in unique(proj$celltype)){
	diseaseMKG(celltype)
}

ff=list.files('DAR_celltype_binomial/','*.bed$')
ff
lapply(ff,function(f) dim(read.table(paste0("DAR_celltype_binomial/",f),sep="\t")))
file <- "./DAR_celltype_binomial/"
ff <- list.files(file,"_DAR.bed$")
peak_lst <- list()
for(i in seq_along(ff)){
	peak_lst[[i]] <- read.table(paste0(file,ff[i]),sep='\t')
    celltype <- gsub('_DAR.bed',"",ff[i])
	peak_lst[[i]]$celltype <- celltype
}
peak_df <- do.call(rbind,peak_lst)
head(peak_df)
colnames(peak_df) <- c("seqnames","start","end","cellname","Log2FC","FDR","celltype")
table(peak_df$celltype)
df <- subset(peak_df,FDR<0.05 & abs(Log2FC)>0.5)

markerpeak <- read.table("./markerPeak_celltype_final/all.celltype.marker.bed",sep='\t',header=T)
head(markerpeak)
table(markerpeak$celltype)
####和每个细胞类型的markerpeak求交集,旨在找出细胞类型特异性的差异motif
df_lst <- list()
k=1
for(i in unique(df$celltype)){
    #peak_gr <- readRDS(paste0("./PeakCalls/",i,"-reproduciblePeaks.gr.rds"))
    peak_gr <- GRanges(subset(markerpeak,celltype==i))
    df_tmp <- GRanges(subset(df,celltype==i))
    x=findOverlaps(df_tmp,peak_gr)
    df_lst[[k]]=as.data.frame(df_tmp[unique(queryHits(x))])
    k=k+1
}
df1 <- do.call(rbind,df_lst)
table(df1$celltype)
head(df1)

##查看每个细胞类型motif富集排行
# find enriched motifs in DARs
#matches <- getMatches(proj, paste(f,"Motif",sep='-'))
matches  <- readRDS("Annotations/Motif-Matches-In-Peaks.rds")
bgdPeaks <- readRDS("Background-Peaks.rds")
bgdPeaks <- SummarizedExperiment::assay(bgdPeaks)
gg_lst <- list()
k=1
for(i in unique(df1$celltype)){
    df_tmp <- subset(df1,celltype==i)
    df_tmp <- GRanges(df_tmp)
    dar <- split(df_tmp,df_tmp$celltype)
    
    enrichList <- lapply(dar, function(d){
      idx <- queryHits(findOverlaps(rowRanges(matches),d))
      ArchR:::.computeEnrichment(matches, idx, c(idx, as.vector(bgdPeaks[idx,])))
    })

    # same wrangling methods as peakAnnoEnrichment
    assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x){
      d <- lapply(seq_along(enrichList), function(y){
        enrichList[[y]][colnames(matches),x,drop=FALSE]
      }) %>% Reduce("cbind",.)
      colnames(d) <- names(enrichList)
      d
    }) %>% SimpleList
    names(assays) <- colnames(enrichList[[1]])
    assays <- rev(assays)
    enrichMotifs <- SummarizedExperiment::SummarizedExperiment(assays=assays)

    df_motif <- data.frame(TF=rownames(enrichMotifs),mlog10Padj=assay(enrichMotifs)[,1])
    df_motif <- df_motif[order(df_motif$mlog10Padj,decreasing=TRUE),]
    df_motif$rank <- seq_len(nrow(df_motif))
    print(head(df_motif))
    df_motif$TF <- gsub('[_].*','',df_motif$TF)
    gg_lst[[k]] <- ggplot(df_motif, aes(rank, mlog10Padj, color = mlog10Padj)) + 
          geom_point(size = 1) +
          ggrepel::geom_label_repel(
                data = df_motif[rev(seq_len(20)), ], aes(x = rank, y = mlog10Padj, label = TF), 
                size = 3,
                nudge_x = 5,
                color = "black"
          ) + theme_ArchR() + 
          ylab("-log10(P-adj) Motif Enrichment") + 
          xlab("Rank Sorted TFs Enriched") +
          scale_color_gradientn(colors = paletteContinuous(set = "comet"))+
          labs(title = paste0(i," DAR enrichMotifs"))
    k=k+1
}
plotPDF(plotList=gg_lst, name = "celltype-DAR-Motifs-Enriched", width = 6, height = 6, ArchRProj = proj, addDOC = FALSE)