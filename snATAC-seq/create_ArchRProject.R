#########################################################################
# File Name: create_ArchRProject.R
# 创建ArchR对象
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
inputdir <- args[1] ###fragments文件所在路径
outdir<-args[2]
setwd(outdir)

genomeAnnotation <- readRDS("../data/Macaca_genomeAnnotation")
geneAnnotation <- readRDS("../data/Macaca_geneAnnotation")

first_name <- list.files(inputdir,"tsv.gz$")
sample_name <- gsub(".fragments.tsv.gz","",first_name)
FragmentFiles <- paste(inputdir,first_name,sep="/")
names(FragmentFiles) <- sample_name

# create ArrowFiles
ArrowFiles <- createArrowFiles(
  inputFiles = FragmentFiles,
  sampleNames = names(FragmentFiles),
  minTSS = 4,   #这个参数不需要过高，后续可以调整
  minFrags = 1000,
  addTileMat = T,
  force= T,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  addGeneScoreMat = TRUE,
  threads=16
)

# add doublet score
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1,
  force=TRUE
)

# create ArchR object
projDM1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "diabetes",
  copyArrows = T,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

# basic QC 
proj_DM_1 <- projDM1
df <- getCellColData(proj_DM_1, select = c("log10(nFrags)", "TSSEnrichment"))
p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = proj_DM_1, addDOC = FALSE)

p1 <- plotGroups(
    ArchRProj = proj_DM_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
p2 <- plotGroups(
    ArchRProj = proj_DM_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p3 <- plotGroups(
    ArchRProj = proj_DM_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "ridges"
   )
p4 <- plotGroups(
    ArchRProj = proj_DM_1, 
    groupBy = "Sample", 
    colorBy = "cellColData", 
    name = "log10(nFrags)",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = proj_DM_1, addDOC = FALSE, width = 10, height = 10)

addArchRThreads(threads = 1)
p1 <- plotFragmentSizes(ArchRProj = proj_DM_1, maxSize=400)
p2 <- plotTSSEnrichment(ArchRProj = proj_DM_1)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = proj_DM_1, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj_DM_1, outputDirectory = "Save-proj_DM", load = FALSE)