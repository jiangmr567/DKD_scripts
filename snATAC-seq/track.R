#########################################################################
# File Name: track.R
# 基因开放性差异
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

###PT and PT_VCAM1的VCAM1基因track图
pt <- proj[proj$celltype %in% c("PT","PT_VCAM1"),]
##motif
tf <- readRDS("Annotations/Motif-Positions-In-Peaks.rds")
head(tf)
tf1 <- tf$NFKB2_199
tf2 <- tf$RELB_312
tf3 <- tf$NFKB1_198
flist <- list()
flist[["peaks"]] <- getPeakSet(pt)
flist[["RELB"]] <- unique(tf2) %>% resize(250, fix="center")
flist[["NFKB1"]] <- unique(tf3) %>% resize(250, fix="center")
flist[["NFKB2"]] <- unique(tf1) %>% resize(250, fix="center")
p <- plotBrowserTrack(
    ArchRProj = pt, 
    groupBy = "celltype",
    geneSymbol = c("VCAM1","TRAF1"),
    #useGroups = atacOrder,
    features = flist,
    #pal =c("#F8766D","#00BFC4"),
    sizes = c(8, 2, 2, 1.5),
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"), 
    loops =getPeak2GeneLinks(
		  ArchRProj = pt, 
		  corCutOff = 0.45, 
		  FDRCutOff = 0.0001,
		  varCutOffATAC = 0.25,
		  varCutOffRNA = 0.25,
		  resolution = 1, 
		  returnLoops = TRUE
	),
    upstream = 20000,
    downstream = 50000,
    tileSize=500,
    minCells=200
)

plotPDF(plotList = p, 
    name = "TRAF1-VCAM1-Tracks.pdf", 
    ArchRProj = pt, 
    addDOC = FALSE, width = 6, height = 6)
    
###PT细胞基因开放性差异分析
PT <- proj[proj$celltype=="PT",]
genes<- c("PCK1","ALDOB","SLC16A9","INSR")

peak <- getPeakSet(PT)
peak_gr <- readRDS("./PeakCalls/PT-reproduciblePeaks.gr.rds")
x=findOverlaps(peak,peak_gr)
peak_filter <-peak[unique(queryHits(x))]
peak_filter

##motif
tf <- readRDS("./Annotations/Motif-Positions-In-Peaks.rds")
tf1 <- tf$HNF4A_582
tf2 <- tf$HNF4G_583
tf3 <- tf$HNF1B_86
flist <- list()
flist[["peaks"]] <- peak_filter
flist[["HNF4A"]] <- unique(tf1) %>% resize(250, fix="center")
flist[["HNF4G"]] <- unique(tf2) %>% resize(250, fix="center")
flist[["HNF1B"]] <- unique(tf3) %>% resize(250, fix="center")

p <- plotBrowserTrack(
    ArchRProj = PT, 
    groupBy = "group",
    geneSymbol = genes,
    #useGroups = atacOrder,
    features = flist,
    sizes = c(10, 2.5, 2.5, 1.5),
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"), 
    loops =getPeak2GeneLinks(
		  ArchRProj = PT, 
		  corCutOff = 0.45, 
		  FDRCutOff = 0.0001,
		  varCutOffATAC = 0.25,
		  varCutOffRNA = 0.25,
		  resolution = 1, 
		  returnLoops = TRUE
	),
    tileSize=500,
    minCells=200
)

plotPDF(plotList = p, 
    name = "PT_maingene-Tracks.pdf", 
    ArchRProj = PT, 
    addDOC = FALSE, width = 6, height = 6)

###Endo细胞基因开放性差异分析 
Endo <- proj[proj$celltype=="Endo",]
genes<- c("TEK","CD74","DOCK4","INSR","NR3C1")
##motif
tf <- readRDS("./Annotations/Motif-Positions-In-Peaks.rds")
tf1 <- tf$ETV1_568
tf2 <- tf$IKZF1_397
tf3 <- tf$GABPA_576
flist <- list()
flist[["peaks"]] <- getPeakSet(Endo)
flist[["ETV1"]] <- unique(tf1) %>% resize(250, fix="center")
flist[["IKZF1"]] <- unique(tf2) %>% resize(250, fix="center")
flist[["GABPA"]] <- unique(tf3) %>% resize(250, fix="center")
p <- plotBrowserTrack(
    ArchRProj = Endo, 
    groupBy = "group",
    geneSymbol = genes,
    #useGroups = atacOrder,
    features = flist,
    sizes = c(10, 2.5, 2.5, 1.5),
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"), 
    loops =getPeak2GeneLinks(
		  ArchRProj = Endo, 
		  corCutOff = 0.45, 
		  FDRCutOff = 0.0001,
		  varCutOffATAC = 0.25,
		  varCutOffRNA = 0.25,
		  resolution = 1, 
		  returnLoops = TRUE
	),
    tileSize=500,
    minCells=200
)
plotPDF(plotList = p, 
    name = "Endo_maingene-Tracks.pdf", 
    ArchRProj = Endo, 
    addDOC = FALSE, width = 6, height = 6)