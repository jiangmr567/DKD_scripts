#########################################################################
# File Name: monocle3.R
# monocle3拟时序分析
#########################################################################
# load library
library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix)
library(irlba)
library(here)
set.seed(1234)

args = commandArgs(T)
outdir<-args[1]
setwd(outdir)

rna <- readRDS("kidney_harmonyfinal.rds")
count_data <- GetAssayData(rna, assay = "RNA",slot = "counts")

gene_metadata <- as.data.frame(rownames(count_data))
colnames(gene_metadata) <- "gene_short_name"
rownames(gene_metadata) <- gene_metadata$gene_short_name
head(gene_metadata)

cds <- new_cell_data_set(as(count_data, "sparseMatrix"),
                         cell_metadata = rna@meta.data,
                         gene_metadata = gene_metadata)
head(cds)
##预处理去批次降维
cds <- preprocess_cds(cds, num_dim = 20)
cds = align_cds(cds, num_dim = 20, alignment_group = "orig.ident")
cds = reduce_dimension(cds,preprocess_method="Aligned")

plot_cells(cds, color_cells_by="celltype", group_label_size = 5)

#subclustering
cds_subset <- choose_cells(cds) 
cds_subset <- cds_subset[,colData(cds_subset)$celltype %in% c("PT","PT_VCAM1")]
plot_cells(cds_subset, color_cells_by = "partition")
cds_subset1 <- choose_cells(cds_subset) 
cds_subset <- cluster_cells(cds_subset1)

cds_subset <- learn_graph(cds_subset)
plot_cells(cds_subset,
           color_cells_by = "celltype",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           label_principal_points = FALSE,
           graph_label_size=5)
pdf("PT_traj_FeaturePlot.pdf",height = 7,width = 7)
plot_cells(cds_subset, genes=c("LRP2","SLC5A2", "SLC7A13","VCAM1"), graph_label_size  = 2,group_label_size = 0)
dev.off()
cds_subset1 <- order_cells(cds_subset)

pdf("PT_traj.pdf",height = 5,width = 5)
plot_cells(cds_subset1,
           color_cells_by = "pseudotime",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           show_trajectory_graph=T)
dev.off()

pdf("PT_plot_genes_in_pseudotime.pdf",height = 2,width = 4.5)
genes <- c("VCAM1","TRAF1","SLC4A4", "SLC34A1","LRP2","SLC13A3")
for (gene in genes) {
  lineage_cds <- cds_subset1[rowData(cds_subset1)$gene_short_name %in% gene,]
  p <- monocle3::plot_genes_in_pseudotime(lineage_cds,
                                     color_cells_by="celltype",
                                     min_expr=1,
                                     cell_size=0.1,
                                     trend_formula = "~ splines::ns(pseudotime, df=10)",
                                     panel_order = gene
  )
  print(p)
}
dev.off()

saveRDS(cds_subset1,"PT_traj_cds.rds")

cds <-cds_subset1
modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
###去除未知基因
genes <- genes[-grep("^ENSMFAG",genes)]    ###73个基因
###画热图
pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
colnames(pt.matrix)<-names(pseudotime(cds)[order(pseudotime(cds))])
pt.matrix[1:3,1:3]
write.csv(pt.matrix,"PT_pseudotime_matrix.csv")

pt.matrix <- read.csv("PT_pseudotime_matrix.csv",row.names=1)

col.anno<-cds@colData[,'celltype'] %>% as.data.frame()
rownames(col.anno) <- rownames(cds@colData)
colnames(col.anno) <- 'celltype'
head(col.anno)

library(pheatmap)
mat <- MinMax(pt.matrix,-2,2)
#pdf("PT_monocle3_heatmap.pdf",width = 6,height = 10)
pheatmap::pheatmap(mat,
                   cluster_cols = F,
                   show_colnames = F,
                   border_color = NA,
                   fontsize_row = 8,
                   annotation_col = col.anno)
#dev.off()

###差异基因交集
deg <- read.csv("kidney_celltype_DM_CTRL_markers.csv",row.names=1)
pt_deg <- subset(deg,celltype %in% c("PT","PT_VCAM1"))
pt_deg_cut <- subset(pt_deg,abs(avg_log2FC)>0.25 & p_val_adj<0.05)
intersect(pt_deg_cut$gene,rownames(pt.matrix))
length(intersect(pt_deg_cut$gene,rownames(pt.matrix)))
length(union(pt_deg_cut$gene,rownames(pt.matrix)))