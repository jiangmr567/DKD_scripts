#########################################################################
# File Name: trajectory_analysis.R
# 使用ArchR进行轨迹分析
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

trajectory <- c("PT", "PT_VCAM1")
proj <- addTrajectory(
    ArchRProj = proj, 
    name = "PT_Traj", 
    groupBy = "celltype",
    trajectory = trajectory, 
    reducedDims = "Harmony",
    embedding = "UMAPHarmony", 
    force = TRUE
)
p <- plotTrajectory(proj, trajectory = "PT_Traj",embedding = "UMAPHarmony", colorBy = "cellColData", name = "PT_Traj")
plotPDF(p, name = "Plot-PT-Traj-UMAP.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
proj <- addImputeWeights(proj)
p1 <- plotTrajectory(proj, trajectory = "PT_Traj",embedding = "UMAPHarmony", colorBy = "GeneScoreMatrix", name = "VCAM1", continuousSet = "horizonExtra")
p2 <- plotTrajectory(proj, trajectory = "PT_Traj",embedding = "UMAPHarmony", colorBy = "GeneScoreMatrix", name = "TRAF1", continuousSet = "horizonExtra")
p3 <- plotTrajectory(proj, trajectory = "PT_Traj",embedding = "UMAPHarmony", colorBy = "GeneScoreMatrix", name = "SLC4A4", continuousSet = "horizonExtra")
p4 <- plotTrajectory(proj, trajectory = "PT_Traj",embedding = "UMAPHarmony", colorBy = "GeneScoreMatrix", name = "LRP2", continuousSet = "horizonExtra")
plotPDF(p1,p2,p3,p4, name = "Plot-PT-maingene-Traj.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

trajMM  <- getTrajectory(ArchRProj = proj, name = "PT_Traj",useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
plotPDF(p1, name = "Plot-PT-motif-Traj.pdf", ArchRProj = proj, addDOC = FALSE, width = 6, height = 10)
trajGSM <- getTrajectory(ArchRProj = proj, name = "PT_Traj", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
tmp <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
colnames(tmp) <- colnames(assay(trajCombined))
assay(trajCombined) <- tmp

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1+ht2, name = "Plot-PT-GSM_MM-Traj.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 10)

trajGIM <- getTrajectory(ArchRProj = proj, name = "PT_Traj", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
corGIM_MM <- correlateTrajectories(trajGIM, trajMM)
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]

trajCombined <- trajGIM2
tmp <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
colnames(tmp) <- colnames(assay(trajCombined))
assay(trajCombined) <- tmp

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
plotPDF(ht1+ht2, name = "Plot-PT-GIM_MM-Traj.pdf", ArchRProj = proj, addDOC = FALSE, width = 10, height = 8)

###足迹分析
motifPositions <- getPositions(proj)
motifs <- c("REL","RELA","RELB","NFKB1","NFKB2")
motifPositions <- getPositions(proj)
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% c("HNF4A.var.2_383","THRB.var.2_460","THRB.var.3_461")]
markerMotifs
seFoot <- getFootprints(
  ArchRProj = proj, 
  positions = motifPositions[markerMotifs], 
  groupBy = "group"
)
seFoot<- seFoot[,grepl("PT_VCAM1", colnames(seFoot))]
plotFootprints(
  seFoot = seFoot,
  ArchRProj = proj, 
  normMethod = "Subtract",
  plotName = "PT_VCAM1-Footprints-Subtract",
  addDOC = FALSE,
  smoothWindow = 5
)