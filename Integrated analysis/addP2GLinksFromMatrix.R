#########################################################################
# File Name: addP2GLinksFromMatrix.R
# 利用scJoint降维矩阵创建peak-gene的link
#########################################################################
# load library
library(Seurat)
library(ArchR)
addArchRThreads(threads = 8)

args = commandArgs(T)
outdir<-args[1]
setwd(outdir)

addPeak2GeneLinksFromMatrix <-  function (ArchRProj, useMatrix = "GeneScoreMatrix",  
    k = 100, knnIteration = 500, overlapCutoff = 0.8, maxDist = 250000, scaleTo = 1e4, 
    log2Norm = TRUE, threads = getArchRThreads(), path='~/d/p2g', file.rna) {
    # path to the Seurat RNA object
    seuratRNA = readRDS(file.rna)
    peakSet <- getPeakSet(ArchRProj)
    geneSet <- ArchR:::.getFeatureDF(getArrowFiles(ArchRProj)[1], useMatrix, threads = threads)
    gs = intersect(geneSet$name, rownames(seuratRNA))
    geneSet <- geneSet[geneSet$name %in% gs,]
    geneStart <- GRanges(geneSet$seqnames, IRanges(geneSet$start, 
        width = 1), name = geneSet$name, idx = geneSet$idx)
    # co-embeddings from scJoint output
    rna  = read.table(file.path(path, 'kidney_exprs_rna_embeddings.txt'))
    atac = read.table(file.path(path, 'kidney_exprs_atac_embeddings.txt'))
    # one cell barcode per line
    rnaNames  = read.csv(file.path(path, 'kidney_rna_cell.csv'))[,2]
    atacNames = read.csv(file.path(path,'kidney_atac_cell.csv'))[,2]
    # select 500 representative points in the co-embedding space
    qry = kmeans(rna, centers=knnIteration, iter.max=10, algorithm="MacQueen")$centers
    knnObj = nabor::knn(data = atac, query = qry, k = k+1)$nn.idx[,-1,drop=FALSE]
    atacObj <- SimpleList(lapply(seq_len(nrow(knnObj)), function(x) atacNames[knnObj[x, ]]))
    knnObj = nabor::knn(data = rna, query = qry, k = k+1)$nn.idx[,-1,drop=FALSE]
    rnaObj <- SimpleList(lapply(seq_len(nrow(knnObj)), function(x) rnaNames[knnObj[x, ]]))
    rm(knnObj,  atac, rna)
    gc()
    peakDF <- mcols(peakSet)
    peakDF$seqnames <- seqnames(peakSet)
    cs  = colnames(seuratRNA)
    mat = Seurat::GetAssayData(seuratRNA, assay = "RNA", slot = "data")
    mat = do.call(cbind,lapply(rnaObj,function(d) Matrix::rowSums(mat[,which(cs %in% d)])))
    groupMatRNA <- mat[geneSet$name,]
    rawMatRNA <- groupMatRNA
    groupMatATAC <- ArchR:::.getGroupMatrix(ArrowFiles = getArrowFiles(ArchRProj), 
        featureDF = peakDF, groupList = atacObj, useMatrix = "PeakMatrix", 
        threads = threads, verbose = FALSE)
    rawMatATAC <- groupMatATAC
    groupMatRNA <- t(t(groupMatRNA)/colSums(groupMatRNA)) * scaleTo
    groupMatATAC <- t(t(groupMatATAC)/colSums(groupMatATAC)) * scaleTo
    if (log2Norm) {
        groupMatRNA <- log2(groupMatRNA + 1)
        groupMatATAC <- log2(groupMatATAC + 1)
    }
    names(geneStart) <- NULL
    seRNA <- SummarizedExperiment(assays = SimpleList(RNA = groupMatRNA, 
        RawRNA = rawMatRNA), rowRanges = geneStart)
    metadata(seRNA)$KNNList <- atacObj
    names(peakSet) <- NULL
    seATAC <- SummarizedExperiment(assays = SimpleList(ATAC = groupMatATAC, 
        RawATAC = rawMatATAC), rowRanges = peakSet)
    metadata(seATAC)$KNNList <- atacObj
    rm(mat, seuratRNA, groupMatRNA, groupMatATAC)
    gc()
    o <- DataFrame(findOverlaps(ArchR:::.suppressAll(resize(seRNA, 2 * maxDist + 1, "center")), 
        resize(rowRanges(seATAC), 1, "center"), ignore.strand = TRUE))
    o$distance <- distance(rowRanges(seRNA)[o[, 1]], rowRanges(seATAC)[o[, 2]])
    colnames(o) <- c("B", "A", "distance")
    o$Correlation <- ArchR:::rowCorCpp(as.integer(o$A), as.integer(o$B), 
        assay(seATAC), assay(seRNA))
    o$VarAssayA <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seATAC)))[o$A]
    o$VarAssayB <- ArchR:::.getQuantiles(matrixStats::rowVars(assay(seRNA)))[o$B]
    o$TStat <- (o$Correlation/sqrt((pmax(1 - o$Correlation^2, 
        1e-17, na.rm = TRUE))/(ncol(seATAC) - 2)))
    o$Pval <- 2 * pt(-abs(o$TStat), ncol(seATAC) - 2)
    o$FDR <- p.adjust(o$Pval, method = "fdr")
    out <- o[, c("A", "B", "Correlation", "FDR", "VarAssayA", "VarAssayB")]
    colnames(out) <- c("idxATAC", "idxRNA", "Correlation", "FDR", "VarQATAC", "VarQRNA")
    mcols(peakSet) <- NULL
    names(peakSet) <- NULL
    metadata(out)$peakSet <- peakSet
    metadata(out)$geneSet <- geneStart
    dir.create("Peak2GeneLinks", showWarnings = FALSE)
    outATAC <- file.path("Peak2GeneLinks", "seATAC-Group-KNN.rds")
    saveRDS(seATAC, outATAC, compress = FALSE)
    outRNA <- file.path("Peak2GeneLinks", "seRNA-Group-KNN.rds")
    saveRDS(seRNA, outRNA, compress = FALSE)
    metadata(out)$seATAC <- outATAC
    metadata(out)$seRNA <- outRNA
    metadata(ArchRProj@peakSet)$Peak2GeneLinks <- out
    saveRDS(ArchRProj, "Save-ArchR-Project.rds")
    ArchRProj
}

proj <- readRDS("Save-ArchR-Project.rds")
path <- "../scJoint/data"
file.rna='../RNA/result/kidney_harmonyfinal.rds'
proj=addPeak2GeneLinksFromMatrix(proj, k = 100, knnIteration = 500, path=path, file.rna=file.rna)