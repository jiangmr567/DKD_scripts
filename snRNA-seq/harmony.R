#########################################################################
# File Name: harmony.R
# 去除批次效应
#########################################################################
# load library
rm(list = ls())
library(Seurat)
library(harmony)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(future)
set.seed(123)  
args = commandArgs(T)
outdir<-args[1]
dimn <- as.numeric(args[2]) #### PC nums

setwd(outdir)
kidney <- readRDS("kidney_preprocessing.rds")
cellinfo <- subset(kidney@meta.data, select = c("orig.ident", "percent.mt", "sample","state"))
kidney <- CreateSeuratObject(kidney@assays$RNA@counts, meta.data = cellinfo)
kidney <- NormalizeData(kidney) %>% FindVariableFeatures(nfeatures = 2000) %>% ScaleData()
kidney <- RunPCA(kidney, npcs=30, verbose=F)

kidney <- RunHarmony(kidney, group.by.vars="orig.ident", max.iter.harmony = 30) 
# group.by.vars参数是设置按哪个分组来整合
# max.iter.harmony设置迭代次数，默认是10。运行RunHarmony结果会提示在迭代多少次后完成了收敛。
# RunHarmony函数中有个lambda参数，默认值是1，决定了Harmony整合的力度。
# lambda值调小，整合力度变大，反之。（只有这个参数影响整合力度，调整范围一般在0.5-2之间）
kidney <- FindNeighbors(kidney, reduction="harmony", dims = 1:dimn)
kidney <- FindClusters(kidney, verbose = FALSE,resolution = c(seq(0.4,1,by=0.1)))
kidney <- RunUMAP(kidney, reduction="harmony", dims = 1:dimn,min.dist = 0.5)
pdf('kidney_harmony_DimPlot.pdf',width = 15,height = 10)
for(i in seq(0.4,1,by=0.1)){
  Idents(object = kidney) <- paste0("RNA_snn_res.",i)
  p <- DimPlot(kidney, reduction = "umap", label = T,label.size = 5)+
    ggtitle(paste0("RNA_snn_res.",i))+theme(plot.title = element_text(hjust = 0.5))
  print(p)
}
DimPlot(kidney, reduction = "umap", group.by = 'sample', split.by = "sample", ncol = 4)
DimPlot(kidney, reduction = "umap", group.by = 'state', split.by = "state", ncol = 2)
dev.off()
kidney@assays
saveRDS(kidney,"kidney_harmony.rds")