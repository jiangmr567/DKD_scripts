#########################################################################
# File Name: celltype_transfer.R
# 利用scJoint方法进行标签转移，并利用降维矩阵整合snRNA-seq数据
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
rna_file <- args[2]  ###RNA对象路径
setwd(outdir)

genomeAnnotation <- readRDS("../data/Macaca_genomeAnnotation")
geneAnnotation <- readRDS("../data/Macaca_geneAnnotation")

proj <- readRDS("Save-ArchR-Project.rds")

###scjoint处理结果
df <- read.csv("../scJoint/data/kidney_umap_embedding.csv",row.names=1)
df <- subset(df,batch=='ATAC')
head(df)
proj$celltype <- df[proj$cellNames,]$celltype

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

###细胞比例百分比图
statistics <- as.data.frame(table(df$state,df$celltype))
colnames(statistics) <- c("state","celltype","cellnum")

cellname <- names(col)
p <- ggplot(data = statistics, mapping = aes(x = state, y = cellnum, 
  fill = factor(celltype,levels=(cellname)))) + 
  geom_bar(stat= 'identity', position = 'fill',width = 0.8)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size = 20))+
  coord_flip()+labs(y="Cell proportion")+
  scale_fill_manual("celltype",values=col)
ggsave("./Plots/atac_state_cellnum_prop.pdf",plot=p,width=10,height=4)


# dimensional reduction
proj <- addIterativeLSI(
    ArchRProj = proj,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 4, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(2), 
        sampleCells = 30000, 
        maxClusters=6,
        n.start = 10
    ), 
    sampleCellsPre = 30000,
    varFeatures = 25000, 
    dimsToUse = 1:25,
    force=T
)

# (Batch correct sample effects)
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
)

# basic clustering 
proj <- addClusters(
    input = proj,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.2,
    force=T
)
set.seed(1)
# UMAP embedding
proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "Harmony", 
    name = "UMAP", 
    nNeighbors = 35, 
    minDist = 0.4, 
    metric = "cosine",force=T
)

###将scJoint得到的降维矩阵及umap坐标替换进ATAC对象
df1 <- df[,c(1,2)]
colnames(df1) <- c("Harmony#UMAP_Dimension_1","Harmony#UMAP_Dimension_2")
proj@embeddings$UMAPHarmony$df <- df1
head(proj@embeddings$UMAPHarmony$df)
###整合RNA数据
emb_df <- read.csv("../scJoint/data/combine_embedding.csv",row.names=1)
df <- emb_df[rownames(proj@reducedDims$Harmony$matDR),]
colnames(df) <- gsub('V','LSI',colnames(df))
proj@reducedDims$Harmony$matDR <- as.matrix(df)
head(proj@reducedDims$Harmony$matDR)

# integration with scRNAseq data
seRNA <- readRDS(paste0(rna_file,"/kidney_harmonyfinal.rds"))

proj <- addGeneIntegrationMatrix(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "Harmony",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    #groupList = groupList,
    groupRNA = "celltype",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)
table(proj$celltype==proj$predictedGroup)
p1 <- plotEmbedding(proj, colorBy = "cellColData", embedding = "UMAPHarmony", name = "predictedGroup")
plotPDF(p1, name = "UMAP_vsRNA.pdf", ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

proj$sample <- gsub("[-_].*","",proj$Sample)
table(proj$sample)
###计算每个样本中位数/平均数
df1 <- aggregate(proj$TSSEnrichment, by=list(type=proj$sample),mean)
df2 <- aggregate(proj$TSSEnrichment, by=list(type=proj$sample),median)
df3 <- aggregate(proj$nFrags, by=list(type=proj$sample),mean)
df4 <- aggregate(proj$nFrags, by=list(type=proj$sample),median)
df5 <- data.frame(table(proj$sample))
df <- data.frame(sample=df1$type,TSS_mean=df1$x,TSS_median=df2$x,nFrags_mean=df3$x,nFrags_median=df4$x,cellnum=df5$Freq)
df$state <- c("CTRL","CTRL","DM","DM")
df

saveRDS(proj,"Save-ArchR-Project.rds")

###参考https://github.com/GreenleafLab/scScalpChromatin/blob/main/Figure_1_Data_Summary.R
source("../utils/plotting_config.R")
source("../utils/misc_helpers.R")
source("../utils/matrix_helpers.R")
source("../utils/archr_helpers.R")

##Ascending LOH ==TAL:thick ascending limb
Asc <- c("SLC12A1","UMOD")
##Proximal tubule cell
PT <- c("SLC34A1","SLC13A3","LRP2","VCAM1")
##Principal cell
PC <- c("FXYD4","GATA2")
##Descending LOH
Des <- c("NLGN1","CLDN1")
##Distal convoluted tubule cell
DCTC <- c("TRPM6","SLC12A3")   
##Connecting tubule cell == CNT
CTC <- c("SLC8A1")
##Endothelial cell
Endo <- c("FLT1","PECAM1")
##Myofibroblast
MyFIB <- c("PRKG1","COL1A2","PDGFRB")
##Type A intercalated cells
ICA <- c("SLC26A7","ATP6V0D2")
##Type B intercalated cells
ICB <- c("SLC26A4")
##Podocyte
Podo <- c("PTPRO","NPHS2")
##Immune
immune <- c("PTPRC","RIPOR2","CD163")

main_markers_gene <- c(PT,Asc,Des,PC,CTC,DCTC,Endo,MyFIB,ICA,ICB,Podo,immune)

GSM_se <- getMatrixFromProject(proj, useMatrix="GeneScoreMatrix")
GSM_mat <- assays(GSM_se)$GeneScoreMatrix
rownames(GSM_mat) <- rowData(GSM_se)$name
seu=CreateSeuratObject(counts=GSM_mat, min.cells=1, min.features=1,meta.data=as.data.frame(GSM_se@colData))
seu=NormalizeData(seu)
saveRDS(seu,'kidney_ArchRtoSeurat.rds')

avgPctMat <- avgAndPctExpressed(GSM_mat[,getCellNames(proj)], proj$celltype, feature_normalize=TRUE, min_pct=0)
avgPctMat <- avgPctMat[avgPctMat$feature %in% main_markers_gene,]
# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 10] <- 0
namedClustAspect <- 1.6
fineClustAspect <- 1.6
grp_order <- rev(levels)
gene_order <- main_markers_gene %>% rev()
p <- dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", 
  xorder=grp_order, yorder=gene_order, cmap=paletteContinuous("solarExtra"), aspectRatio=namedClustAspect)
pdf("./Plots/atac_celltype_marker_dotplot.pdf",width=6,height=10)
p
dev.off()
###利用转成seurat对象画UMAP
library(Seurat)
kidney <- readRDS("kidney_ArchRtoSeurat.rds")
### 标准化数据
kidney <- NormalizeData(kidney) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
### PCA
kidney <- RunPCA(kidney,verbose=F)
#kidney <- FindNeighbors(kidney,dims = 1:25)
#kidney <- FindClusters(kidney, verbose = FALSE,resolution =0.3)
kidney <- RunUMAP(kidney, dims = 1:25)

###scjoint处理结果
df <- read.csv("../scJoint/data/kidney_umap_embedding.csv",row.names=1)
df <- subset(df,batch=='ATAC')
df2 <- df[,c(1,2)]
colnames(df2) <-c("UMAP_1","UMAP_2")
kidney@reductions$umap@cell.embeddings <- as.matrix(df2)
Idents(kidney) <- "celltype"
levels(kidney) <- rev(levels)
levels(kidney)
kidney$celltype <- factor(kidney$celltype,levels=rev(levels))
saveRDS(kidney,"kidney_ArchRtoSeurat.rds")
pdf('./Plots/ATAC_celltype_DimPlot.pdf',width = 8,height = 6)
DimPlot(kidney,cols=col,reduction = "umap", label = T,group.by = 'celltype',raster=FALSE)
dev.off()

pdf('./Plots/ATAC_DM_CTRL.pdf',width = 13,height = 6)
DimPlot(kidney,cols=col,reduction = "umap", label = T,group.by = 'celltype',split.by="state",raster=FALSE)
dev.off()

###聚类基因表达图
q<-FeaturePlot(kidney,features = main_markers_gene,combine = F,cols =paletteContinuous("solarExtra"),raster=FALSE,order=T)
for(i in 1:length(main_markers_gene)){
  q[[i]]<-q[[i]]+xlab("")+ylab("")+NoLegend()
}
pdf("./Plots/kidney_genescore_FeaturePlot.pdf", onefile = TRUE,width = 10,height = 10)  
m<-1
for(i in 1:(floor(length(main_markers_gene)/16)+1)){
  p<-plot_grid(plotlist=q[m:(i*16)])
  m<-i*16+1
  print(p)
}
dev.off()

###ATAC质控
atac_ccd <- proj@cellColData %>% as.data.frame()
dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)
atac_col <- colors
names(atac_col) <- unique(proj$Sample)
atac_col
names(atac_col)
# ATAC TSS / cell
p <- (
    ggplot(atac_ccd, aes(x=Sample, y=TSSEnrichment, fill=Sample))
    + geom_violin(aes(fill=Sample), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
    + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
    + scale_color_manual(values=atac_col, limits=names(atac_col), name="Sample", na.value="grey")
    + scale_fill_manual(values=atac_col)
    + guides(fill=guide_legend(title="Sample"), 
      colour=guide_legend(override.aes = list(size=5)))
    + xlab("")
    + ylab("TSS Enrichment")
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    + scale_y_continuous(limits=c(0,25), expand = c(0, 0))
)

# ATAC log10 nFrags / cell
p1 <- (
    ggplot(atac_ccd, aes(x=Sample, y=log10(nFrags), fill=Sample))
    + geom_violin(aes(fill=Sample), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
    + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
    + scale_color_manual(values=atac_col, limits=names(atac_col), name="Sample", na.value="grey")
    + scale_fill_manual(values=atac_col)
    + guides(fill=guide_legend(title="Sample"), 
      colour=guide_legend(override.aes = list(size=5)))
    + xlab("")
    + ylab("log10 nFrags")
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    + scale_y_continuous(limits=c(0,5), expand = c(0, 0))
)
pdf("./Plots/ATAC_qc_violin.pdf", width=10, height=4)
p
p1
dev.off()

##marker基因track图
markers <- main_markers_gene
uniq_symbol <- unique(proj@geneAnnotation$genes$symbol)
Genes <- intersect(markers, uniq_symbol)
df <- geneAnnotation$gene[geneAnnotation$gene$symbol %in% Genes,]
df1 <- as.data.frame(df)
df1 <- GRanges(df1)
plotRegions <- resize(df1, width=width(df1) + 0.05*width(df1), fix="center")
###按照提供基因顺序画
plotRegions$symbol <- factor(plotRegions$symbol,levels=markers)
plotRegions <- plotRegions[order(plotRegions$symbol),]
levels=c("Immune", "PODO", "ICB", "ICA", "Myofib", "Endo", "DCT", 
         "CNT", "PC", "PEC", "TAL","PT_VCAM1","PT")
celltype_order <- rev(levels)
# Tracks of genes (scalp):
p <- plotBrowserTrack(
    ArchRProj = proj, 
    groupBy = "celltype", 
    useGroups = celltype_order,
    pal = col,
    sizes = c(12, 3, 4, 1),
    region=plotRegions,
    geneSymbol = Genes, 
    loops = NULL,
    features =NULL,
    tileSize=500,
    minCells=200
)
plotPDF(plotList = p, 
    name = "markergenes-Tracks.pdf", 
    ArchRProj = proj, 
    addDOC = FALSE, width = 6, height = 8)