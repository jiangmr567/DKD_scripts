#########################################################################
# File Name: data_to_h5.R
# 将基因表达矩阵和genescore矩阵转换成h5文件
#########################################################################
library(ArchR)
library(parallel)
library(Seurat)
library(cowplot)
library(harmony)
library(ggplot2)
library(tibble)
library(dplyr)
library(Matrix)
library(rhdf5)
library(HDF5Array)

args = commandArgs(T)
outdir<-args[1]
input_rna <- args[2] ###snRNA-seq数据路径
input_atac <- args[3] ###snATAC-seq数据路径
setwd(outdir)

#' This function will generate h5 files for a list expression matrices as input of process_db.py file
write_h5_scJoint <- function(exprs_list, h5file_list) {
  for (i in seq_along(exprs_list)) {
    if (file.exists(h5file_list[i])) {
      warning("h5file exists! will rewrite it.")
      system(paste("rm", h5file_list[i]))
    }
    
    h5createFile(h5file_list[i])
    h5createGroup(h5file_list[i], "matrix")
    writeHDF5Array(t((exprs_list[[i]])), h5file_list[i], name = "matrix/data")
    h5write(rownames(exprs_list[[i]]), h5file_list[i], name = "matrix/features")
    h5write(colnames(exprs_list[[i]]), h5file_list[i], name = "matrix/barcodes")
    print(h5ls(h5file_list[i]))
    
  }
}

#' This function will generate csv files for a list of cell types as input of process_db.py file
write_csv_scJoint <- function(cellType_list, csv_list) {
  
  for (i in seq_along(cellType_list)) {
    
    if (file.exists(csv_list[i])) {
      warning("csv_list exists! will rewrite it.")
      system(paste("rm", csv_list[i]))
    }
    
    names(cellType_list[[i]]) <- NULL
    write.csv(cellType_list[[i]], file = csv_list[i])
    
  }
}

rna <- readRDS(paste0(input_rna,"/kidney_harmonyfinal.rds"))
atac <- loadArchRProject(paste0(input_atac,"/Save-proj_DM"))
write.csv(colnames(rna@assays$RNA@counts),"kidney_rna_cell.csv")

mtx=getMatrixFromProject(atac, useMatrix="GeneScoreMatrix")
count=mtx@assays@data[[1]]
rownames(count)=rowData(mtx)$name
print(count[1:3,1:3])
print(table(colnames(count)==atac$cellNames))
write.csv(colnames(count),"kidney_atac_cell.csv")

rna_gene <- rownames(rna@assays$RNA@counts)
atac_gene <- rownames(count)
markers <- read.csv(paste0(input_rna,"/kidney_celltype_markers.csv"))
markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC) -> top50_markers
gene <- Reduce(intersect,list(rna_gene,atac_gene,rna@assays$RNA@var.features))
gene <- union(intersect(top50_markers$gene,atac_gene),gene)
print(length(gene))
expres_rna <- Matrix(rna@assays$RNA@counts[gene,],sparse=F)
expres_atac <- Matrix(count[gene,],sparse=F)
write_h5_scJoint(exprs_list = list(rna = expres_rna,
								   atac = expres_atac), 
				 h5file_list = c("kidney_exprs_rna.h5", 
								 "kidney_exprs_atac.h5"))
write_csv_scJoint(cellType_list =  list(rna = rna$celltype),
				 csv_list = c("kidney_celltype_rna.csv"))

