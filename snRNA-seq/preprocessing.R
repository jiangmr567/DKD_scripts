#########################################################################
# File Name: preprocessing.R
# 数据预处理，包括质控和去双胞
#########################################################################
# load library
rm(list = ls())
library(Seurat)
library(tidyverse)
library(tidyr)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(future)
library(DoubletFinder)
set.seed(123) 
args = commandArgs(T)
indir <- args[1]        #### soupX_matrix
outdir<-args[2]
dimn <- as.numeric(args[3]) #### PC nums
res <- as.numeric(args[4]) #### resolution

setwd(outdir)
##==合并数据集==##
lst <- list()
i=1
# 以下代码会把每个样本的数据创建一个seurat对象，并存放到列表lst里
for(file in list.files(indir)){
  try({expression_matrix <- Read10X(paste0(indir,'/',file))
    lst[[i]] <- CreateSeuratObject(counts = expression_matrix, project = file, min.cells = 3, min.features = 200)
    # 给细胞barcode加个前缀，防止合并后barcode重名
    lst[[i]] <- RenameCells(lst[[i]], add.cell.id = file)  
    # 计算线粒体基因比例
    mt.gene <- c("ATP6","ATP8","COX1","COX2","COX3","CYTB","ND1","ND2","ND3","ND4","ND4L","ND5","ND6")
    mt.gene <- intersect(mt.gene,rownames(lst[[i]]))
    lst[[i]][["percent.mt"]] <- PercentageFeatureSet(lst[[i]], features = mt.gene)
    i=i+1
  })
}
kidney_initial <-  merge(x=lst[[1]],y=lst_qc[2:length(lst)])
saveRDS(kidney_initial,"kidney_initial.rds")

## 质控，去双胞
Find_doublet <- function(data){
  sweep.res.list <- paramSweep_v3(data, PCs = 1:10, sct = FALSE)
  # 使用log标准化，sct参数设置为 sct = F（默认 ）,如使用SCT标准化方法，设置为T
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)   # 可以看到最佳参数的点
  DoubletRate = ncol(data)*8*1e-6   # 按每增加1000个细胞，双细胞比率增加千分之8来计算
  nExp_poi <- round(DoubletRate*ncol(data))
  p <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric() # 提取最佳pk值
  data <- doubletFinder_v3(data, PCs = 1:10, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
  data
}

doublet=list()
lst_qc=list() 
for (i in 1:length(x = lst)) {
  lst[[i]] <- subset(lst[[i]], subset = nFeature_RNA < 5000 & nFeature_RNA > 500 & nCount_RNA < 20000 & percent.mt < 10)
  lst_qc[[i]] <- lst[[i]]
  lst[[i]] <- NormalizeData(object = lst[[i]], verbose = FALSE)
  lst[[i]] <- FindVariableFeatures(object = lst[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  lst[[i]] <- ScaleData(lst[[i]])
  lst[[i]] <- RunPCA(lst[[i]],npcs=30)
  lst[[i]] <- RunUMAP(lst[[i]], dims = 1:10)
  lst[[i]] <- Find_doublet(lst[[i]])
  doublet[[i]] <- lst[[i]]@meta.data[,c(-2,-3,-4,-5)]
  lst[[i]] <- subset(lst[[i]],subset=doublet_info=="Singlet")
}
# 保存双胞信息
doublet_df=do.call(rbind,doublet)
write.table(doublet_df,"Doublet_info_kidney.txt",sep="\t",quote=FALSE)

# qc后但未去双胞
kidney_qc <-  merge(x=lst_qc[[1]],y=lst_qc[2:length(lst_qc)])
saveRDS(kidney_qc,"kidney_qc.rds")

# 去双胞后
kidney <-  merge(x=lst[[1]],y=lst[2:length(lst)])
kidney@meta.data$sample <- substring(kidney@meta.data$orig.ident,1,2)
kidney@meta.data$state <- "1"
kidney@meta.data$state[grep(pattern="^M[1,2]",kidney@meta.data$orig.ident)] <- "CTRL"
kidney@meta.data$state[grep(pattern="^M[7,9]",kidney@meta.data$orig.ident)] <- "DM"
dim(kidney)
table(kidney@meta.data$sample)
table(kidney@meta.data$state)
table(kidney@meta.data$orig.ident)

pdf('kidney_VlnPlot_after_preprocessing.pdf',width = 10,height = 8)
VlnPlot(kidney,features = c("nFeature_RNA"),pt.size=0) + NoLegend()
VlnPlot(kidney,features = c("nCount_RNA"),pt.size=0) + NoLegend()
VlnPlot(kidney,features = c("percent.mt"),pt.size=0) + NoLegend()
dev.off()

saveRDS(kidney,"kidney_preprocessing.rds")
