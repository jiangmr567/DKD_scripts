#########################################################################
# File Name: celltype_annotation.R
# 细胞类型注释
#########################################################################
# load library
rm(list = ls())
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
set.seed(123) 

args = commandArgs(T)
outdir<-args[1]
setwd(outdir)

kidney <- readRDS("kidney_harmony.rds")
###计算每个样本中位数/平均数
df1 <- aggregate(kidney1$nCount_RNA, by=list(type=kidney1$sample),mean)
df2 <- aggregate(kidney1$nCount_RNA, by=list(type=kidney1$sample),median)
df3 <- aggregate(kidney1$nFeature_RNA, by=list(type=kidney1$sample),mean)
df4 <- aggregate(kidney1$nFeature_RNA, by=list(type=kidney1$sample),median)
df5 <- data.frame(table(kidney1$sample))
df <- data.frame(sample=df1$type,UMI_mean=df1$x,UMI_median=df2$x,ngene_mean=df3$x,ngene_median=df4$x,cellnum=df5$Freq)
df
Idents(object = kidney) <- "RNA_snn_res.0.6"

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
##low quality
low <- c("EYS")

main_markers_gene <- c(PT,Asc,Des,PC,CTC,DCTC,Endo,MyFIB,ICA,ICB,Podo,immune,low)
#聚类基因表达图
DefaultAssay(object = kidney) <- "RNA"
q<-FeaturePlot(kidney,features = main_markers_gene,combine = F,cols = c("grey","red"))
for(i in 1:length(main_markers_gene)){
  q[[i]]<-q[[i]]+xlab("")+ylab("")+NoLegend()
}
pdf("kidney_harmonyFeaturePlot.pdf", onefile = TRUE,width = 10,height = 10)  
m<-1
for(i in 1:(floor(length(main_markers_gene)/16)+1)){
  p<-plot_grid(plotlist=q[m:(i*16)])
  m<-i*16+1
  print(p)
}
dev.off()

#基因表达小提琴图
q<-VlnPlot(kidney, features =main_markers_gene, pt.size = 0, combine = FALSE)
for(i in 1:length(main_markers_gene)){
  q[[i]]<-q[[i]]+theme(axis.text.x = element_text(size = 8))+xlab("")+ylab("")+NoLegend()
}
pdf("kidney_harmonyVlnPlot.pdf", onefile = TRUE,width = 12,height = 9)
m<-1
for(i in 1:(floor(length(main_markers_gene)/9)+1)){
  p<-plot_grid(plotlist=q[m:(i*9)])
  m<-i*9+1
  print(p)
}
dev.off()

DefaultAssay(object = kidney) <- "RNA"
## 细胞类型注释
new.cluster.ids <- c("PT","PC","TAL","PT","PT","DCT","TAL","PT","Endo","low quality",
					 "CNT", "ICA", "Myofib","PEC","PT","PT_VCAM1","Endo","DCT",
					 "Immune","Myofib","PODO","ICB","TAL","PC","Immune","Endo")
names(new.cluster.ids) <- levels(kidney)
new.cluster.ids
# 更改细胞聚类名字
kidney <- RenameIdents(kidney, new.cluster.ids)
kidney@meta.data$celltype <- Idents(kidney)
levels(kidney)
pdf('kidney_harmony_annotation.pdf',width = 15,height = 10)
DimPlot(kidney, reduction = "umap", label = T)
dev.off()

# 保存结果
saveRDS(kidney,"kidney_harmonyfile.rds")

#去除低质量细胞
DefaultAssay(object = kidney) <- "RNA"
dim(kidney)
kidney <- subset(kidney,subset=celltype!="low quality")
saveRDS(kidney,"kidney_harmonyfinal.rds")

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
levels(kidney) <- rev(levels)
levels(kidney)
kidney$celltype <- factor(kidney$celltype,levels=rev(levels))
unique(kidney$celltype)

###细胞类型UMAP
pdf('kidney_final_DimPlot.pdf',width = 8,height = 6)
DimPlot(kidney,cols=col,reduction = "umap", label = T,group.by = 'celltype',pt.size=0.3)
dev.off()

pdf('kidney_state_DimPlot.pdf',width = 13,height = 6)
DimPlot(kidney,cols=col,reduction = "umap", label = T,group.by = 'celltype',split.by="state",pt.size=0.3)
dev.off()

#去除低质量细胞后计算每个样本中位数/平均数
df1 <- aggregate(kidney$nCount_RNA, by=list(type=kidney$sample),mean)
df2 <- aggregate(kidney$nCount_RNA, by=list(type=kidney$sample),median)
df3 <- aggregate(kidney$nFeature_RNA, by=list(type=kidney$sample),mean)
df4 <- aggregate(kidney$nFeature_RNA, by=list(type=kidney$sample),median)
df5 <- data.frame(table(kidney$sample))
df <- data.frame(sample=df1$type,UMI_mean=df1$x,UMI_median=df2$x,ngene_mean=df3$x,ngene_median=df4$x,cellnum=df5$Freq)
df$state <- c("CTRL","CTRL","DM","DM")
df

###细胞比例百分比图
statistics <- as.data.frame(table(kidney$state,kidney$celltype))
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
ggsave("state_cellnum_prop.pdf",plot=p,width=10,height=4)

###参考https://github.com/GreenleafLab/scScalpChromatin/blob/main/Figure_1_Data_Summary.R
source("../utils/plotting_config.R")
source("../utils/misc_helpers.R")
source("../utils/matrix_helpers.R")
source("../utils/archr_helpers.R")
main_markers_gene <- c(PT,Asc,Des,PC,CTC,DCTC,Endo,MyFIB,ICA,ICB,Podo,immune)
###聚类基因表达图
q<-FeaturePlot(kidney,features = main_markers_gene,combine = F,cols =paletteContinuous("blueYellow"),order=T)
for(i in 1:length(main_markers_gene)){
  q[[i]]<-q[[i]]+xlab("")+ylab("")+NoLegend()
}
pdf("kidney_harmonyFeaturePlot.pdf", onefile = TRUE,width = 10,height = 10)  
m<-1
for(i in 1:(floor(length(main_markers_gene)/16)+1)){
  p<-plot_grid(plotlist=q[m:(i*16)])
  m<-i*16+1
  print(p)
}
dev.off()

###marker基因点图
count_mat <- GetAssayData(object=kidney, slot="counts")
avgPctMat <- avgAndPctExpressed(count_mat, kidney$celltype, feature_normalize=TRUE, min_pct=0)
head(avgPctMat)
avgPctMat <- avgPctMat[avgPctMat$feature %in% main_markers_gene,]
# Threshold min pct
avgPctMat$pctExpr[avgPctMat$pctExpr < 5] <- 0
namedClustAspect <- 1.6
fineClustAspect <- 1.6
grp_order <- rev(levels)
gene_order <- main_markers_gene %>% rev()
p <- dotPlot(avgPctMat, xcol="grp", ycol="feature", color_col="avgExpr", size_col="pctExpr", 
  xorder=grp_order, yorder=gene_order, cmap=paletteContinuous("blueYellow"), aspectRatio=namedClustAspect)
pdf("celltype_marker_dotplot.pdf",width=6,height=10)
p
dev.off()

rna_ccd <- kidney@meta.data
rna_ccd$Sample <- rna_ccd$orig.ident
dodge_width <- 0.75
dodge <- position_dodge(width=dodge_width)
rna_col <- colors
names(rna_col) <- unique(kidney$orig.ident)
rna_col
names(rna_col)
# RNA nUMIs / cell
p <- (
    ggplot(rna_ccd, aes(x=Sample, y=nCount_RNA/1000, fill=Sample))
    + geom_violin(aes(fill=Sample), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
    + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
    + scale_color_manual(values=rna_col, limits=names(rna_col), name="Sample", na.value="grey")
    + scale_fill_manual(values=rna_col)
    + guides(fill=guide_legend(title="Sample"), 
      colour=guide_legend(override.aes = list(size=5)))
    + xlab("")
    + ylab("Number of UMIs(x1000)")
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    + scale_y_continuous(limits=c(0,15), expand = c(0, 0))
)
# RNA nGenes / cell
p1 <- (
    ggplot(rna_ccd, aes(x=Sample, y=nFeature_RNA/1000, fill=Sample))
    + geom_violin(aes(fill=Sample), alpha=0.5, adjust = 1.0, scale='width', position=dodge)
    + geom_boxplot(alpha=0.5, width=0.25, outlier.shape = NA)
    + scale_color_manual(values=rna_col, limits=names(rna_col), name="Sample", na.value="grey")
    + scale_fill_manual(values=rna_col)
    + guides(fill=guide_legend(title="Sample"), 
      colour=guide_legend(override.aes = list(size=5)))
    + xlab("")
    + ylab("Number of genes(x1000)")
    + theme_BOR(border=FALSE)
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            #aspect.ratio = aspectRatio, # What is the best aspect ratio for this chart?
            legend.position = "none", # Remove legend
            axis.text.x = element_text(angle = 90, hjust = 1))
    + scale_y_continuous(limits=c(0,4), expand = c(0, 0))
)
pdf("RNA_qc_violin.pdf", width=10, height=4)
p
p1
dev.off()