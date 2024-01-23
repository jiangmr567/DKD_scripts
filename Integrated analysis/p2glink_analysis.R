#########################################################################
# File Name: p2glink_analysis.R
# peak-gene link分析
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
###参考https://github.com/GreenleafLab/scScalpChromatin/blob/main/Figure_2_Linked_Peaks.R
source("../utils/plotting_config.R")
source("../utils/misc_helpers.R")
source("../utils/matrix_helpers.R")
source("../utils/archr_helpers.R")

plotDir <- "p2gLink_plots"
dir.create(plotDir, showWarnings = FALSE, recursive = TRUE)

colors = c(
  "#1B9E77", "#80B1D3", "#BEAED4", "#FB8072", "#386CB0", "#F0027F",
  "#E78AC3", "#BF5B17", "#7FC97F", "#FF7F00", "#7570B3", "#E6AB02", 
  "#FF0000", "#A6761D", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#E5D8BD", "#CAB2D6", "#6A3D9A", "#B15928",
  "#FBB4AE", "#B3CDE3", "#BC80BD", "#CCEBC5", "#DECBE4", "#FED9A6",
  "#FFFFB3", "#FDDAEC", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4",
  "#E6F5C9", "#FFF2AE", "#F1E2CC", "#E41A1C", "#377EB8", "#984EA3",
  "#D95F02", "#F781BF", "#FFED6F", "#66C2A5", "#FC8D62", "#8DA0CB", 
  "#A6D854", "#FFD92F", "#E5C494", "#8DD3C7", "#A65628", "#FCCDE5",
  "#BEBADA", "#FDC086", "#A6CEE3", "#FDB462")
levels=c("Immune", "PODO", "ICB", "ICA", "Myofib", "Endo", "DCT", 
         "CNT", "PC", "PEC", "TAL","PT_VCAM1","PT")
names(colors) <- rev(levels)
col <- colors[1:13]

# Get all peaks
allPeaksGR <- getPeakSet(proj)
allPeaksGR$peakName <- (allPeaksGR %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})
names(allPeaksGR) <- allPeaksGR$peakName
head(allPeaksGR)

# Prepare lists to store peaks, p2g links, loops, coaccessibility
plot_loop_list <- list()
plot_loop_list[["kidney"]] <- getPeak2GeneLinks(proj, resolution = 100)[[1]]
peak2gene_list <- list()
p2gGR <- getP2G_GR(proj, corrCutoff=NULL, varCutoffATAC=-Inf, varCutoffRNA=-Inf, filtNA=FALSE)
p2gGR$source <- "kidney"
peak2gene_list[["kidney"]] <- p2gGR
full_p2gGR <- as(peak2gene_list, "GRangesList") %>% unlist()
saveRDS(full_p2gGR, file="p2gLink_plots/multilevel_p2gGR.rds") # NOT merged or correlation filtered
saveRDS(plot_loop_list, file="p2gLink_plots/multilevel_plot_loops.rds")

##########################################################################################
# 过滤冗余峰与基因的链接
##########################################################################################

# Get metadata from full project to keep for new p2g links
originalP2GLinks <- metadata(proj@peakSet)$Peak2GeneLinks
p2gMeta <- metadata(originalP2GLinks)

# 折叠多余的 p2gLinks：
full_p2gGR <- full_p2gGR[order(full_p2gGR$Correlation, decreasing=TRUE)]
filt_p2gGR <- full_p2gGR[!duplicated(paste0(full_p2gGR$peakName, "_", full_p2gGR$symbol))] %>% sort()

# Reassign full p2gGR to archr project
new_p2g_DF <- mcols(filt_p2gGR)[,c(1:6)]
metadata(new_p2g_DF) <- p2gMeta
metadata(proj@peakSet)$Peak2GeneLinks <- new_p2g_DF

##########################################################################################
# Plot Peak2Gene heatmap
##########################################################################################
nclust <- 25 
p <- plotPeak2GeneHeatmap(
  proj, 
  corCutOff = 0.45, 
  groupBy="celltype", 
  nPlot = 1000000, returnMatrices=FALSE, 
  k=nclust, seed=1, palGroup=col
)
pdf(paste0(plotDir, "/kidney_peakToGeneHeatmap.pdf"), width=16, height=10)
p
dev.off()

# Need to force it to plot all peaks if you want to match the labeling when you 'returnMatrices'.
p2gMat <- plotPeak2GeneHeatmap(
  proj, 
  corCutOff = 0.45, 
  groupBy="celltype",
  nPlot = 1000000, returnMatrices=TRUE, 
  k=nclust, seed=1)
# Get association of peaks to clusters
kclust_df <- data.frame(
  kclust=p2gMat$ATAC$kmeansId,
  peakName=p2gMat$Peak2GeneLinks$peak,
  gene=p2gMat$Peak2GeneLinks$gene
  )
  
# Fix peakname
kclust_df$peakName <- sapply(kclust_df$peakName, function(x) strsplit(x, ":|-")[[1]] %>% paste(.,collapse="_"))
head(kclust_df)
# Get motif matches
matches <- getMatches(proj, "Motif")
r1 <- SummarizedExperiment::rowRanges(matches)
rownames(matches) <- paste(seqnames(r1),start(r1),end(r1),sep="_")
matches <- matches[names(allPeaksGR)]
clusters <- unique(kclust_df$kclust) %>% sort()
clusters
enrichList <- lapply(clusters, function(x){
  cPeaks <- kclust_df[kclust_df$kclust == x,]$peakName %>% unique()
  ArchR:::.computeEnrichment(matches, which(names(allPeaksGR) %in% cPeaks), seq_len(nrow(matches)))
  }) %>% SimpleList
names(enrichList) <- clusters

# Format output to match ArchR's enrichment output
assays <- lapply(seq_len(ncol(enrichList[[1]])), function(x){
    d <- lapply(seq_along(enrichList), function(y){
        enrichList[[y]][colnames(matches),x,drop=FALSE]
      }) %>% Reduce("cbind",.)
    colnames(d) <- names(enrichList)
    d
  }) %>% SimpleList
names(assays) <- colnames(enrichList[[1]])
assays <- rev(assays)
res <- SummarizedExperiment::SummarizedExperiment(assays=assays)

formatEnrichMat <- function(mat, topN, minSig, clustCols=TRUE){
  plotFactors <- lapply(colnames(mat), function(x){
    ord <- mat[order(mat[,x], decreasing=TRUE),]
    ord <- ord[ord[,x]>minSig,]
    rownames(head(ord, n=topN))
  }) %>% do.call(c,.) %>% unique()
  pMat <- mat[plotFactors,]
  prettyOrderMat(pMat, clusterCols=clustCols)$mat
}

pMat <- formatEnrichMat(assays(res)$mlog10Padj, 5, 10, clustCols=FALSE)
# Save maximum enrichment
tfs <- strsplit(rownames(pMat), "_") %>% sapply(., `[`, 1)
rownames(pMat) <- paste0(tfs, " (", apply(pMat, 1, function(x) floor(max(x))), ")")

pMat <- apply(pMat, 1, function(x) x/max(x)) %>% t()
pdf(paste0(plotDir, "/kidney_enrichedMotifs_kclust.pdf"), width=12, height=12)
ht_opt$simple_anno_size <- unit(0.25, "cm")
hm <- BORHeatmap(
  pMat, 
  limits=c(0,1), 
  clusterCols=FALSE, clusterRows=FALSE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = cmaps_BOR$comet,
  #top_annotation = ta,
  row_names_side = "left",
  width = ncol(pMat)*unit(0.5, "cm"),
  height = nrow(pMat)*unit(0.33, "cm"),
  border_gp=gpar(col="black"), # Add a black border to entire heatmap
  legendTitle="Norm.Enrichment -log10(P-adj)[0-Max]"
  )
draw(hm)
dev.off()

# GO enrichments of top N genes per cluster 
# ("Top" genes are defined as having the most peak-to-gene links)
source("../utils/GO_wrappers.R")
kclust <- unique(kclust_df$kclust) %>% sort()
all_genes <- kclust_df$gene %>% unique() %>% sort()

# Save table of top linked genes per kclust
nGOgenes <- 200
topKclustGenes <- lapply(kclust, function(k){
  kclust_df[kclust_df$kclust == k,]$gene %>% getFreqs() %>% head(nGOgenes) %>% names()
  }) %>% do.call(cbind,.)
write.table(topKclustGenes, paste0(plotDir, "/top200_genes_kclust.tsv"), quote=FALSE, sep='\t', row.names = FALSE, col.names=TRUE)

GOresults <- lapply(kclust, function(k){
  message(sprintf("Running GO enrichments on k cluster %s...", k))
  clust_genes <- topKclustGenes[,k]
  upGO <- rbind(
    calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="BP") 
    #calcTopGo(all_genes, interestingGenes=clust_genes, nodeSize=5, ontology="MF")
    #calcTopGo(all_genes, interestingGenes=upGenes, nodeSize=5, ontology="CC")
    )
  upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]
  })

names(GOresults) <- paste0("cluster_", kclust)

# Plots of GO term enrichments:
pdf(paste0(plotDir, "/kclust_GO_3termsBPonlyBarLim.pdf"), width=10, height=2.5)
for(name in names(GOresults)){
    goRes <- GOresults[[name]]
    if(nrow(goRes)>1){
      print(topGObarPlot(goRes, cmap = cmaps_BOR$comet, 
        nterms=3, border_color="black", 
        barwidth=0.85, title=name, barLimits=c(0, 15)))
    }
}
dev.off()