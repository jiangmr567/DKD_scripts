#########################################################################
# File Name: GWAS_plot.R
# GWAS富集可视化
#########################################################################
# load library
library(dplyr)

args = commandArgs(T)
outdir<-args[1]
sumstat_path <- args[2]
setwd(outdir)

celltype=rev(c("Immune", "PODO", "ICB", "ICA", "Myofib", "Endo", "DCT", "CNT","PEC", "PC", "TAL", "PT_VCAM1","PT"))
sumstats_df <- read.table(paste0(sumstat_path,"/diabets_macaca.sumstats.file.list.txt"),sep='\t')
all.sumstats.name = sort(sumstats_df$V1)

#生成每种细胞类型的p值文件
ff <- list.files("../",".txt$")
lst <- list()
i=1
for(f in ff){
	sumstat.name <- gsub(".regressions.*","",f)
	lst[[i]] <- read.table(paste0("../",f),sep = "\t",header=T)
	lst[[i]]$sumstat.names = sumstat.name
	i=i+1
}

#将所有细胞类型的p值文件合并
all.pvalue <- do.call(rbind,lst)
write.csv(all.pvalue,"all_celltype_sumstat.csv")

all.pvalue <- read.csv("all_celltype_sumstat.csv")
all.pvalue.s = all.pvalue[,-c(1,3,4)]
colnames(all.pvalue.s) = c("celltype","p.value","sumstats")
all.pvalue.s$log.p = -log(all.pvalue.s$p.value)
all.pvalue.s = all.pvalue.s[,-2]
#长表转成宽表格
library(tidyr)
all.pvalue.s.kuan = spread(all.pvalue.s,sumstats,log.p)
rownames(all.pvalue.s.kuan) =all.pvalue.s.kuan$celltype
all.pvalue.s.kuan= all.pvalue.s.kuan[,-1]

library(Seurat)
matrix=MinMax(data = all.pvalue.s.kuan, min = 0, max = 4)

###画celltype和sumstats相关性热图
library(pheatmap)
library(ggplot2)

pheatmap=pheatmap(t(matrix),scale='none',cluster_row = TRUE,cluster_cols=TRUE)
ggsave(plot=pheatmap,"peak.sumstats.pvalue.pheatmap.pdf",width=8,height=5)


##GWAS富集柱状图
library(stringr)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tibble)
library(ComplexHeatmap)
library(Seurat)

df <- read.csv("all_celltype_sumstat.csv",row.names=1)
colnames(df) <- c('celltype','Coefficient','Coefficient_std_error','Coefficient_P_value','sumstats')
celltype=rev(c("Immune", "PODO", "ICB", "ICA", "Myofib", "Endo", "DCT", "CNT","PEC", "PC", "TAL", "PT_VCAM1","PT"))
sumstats <- c("chronic_kidney_disease","Microalbuminuria","MW_eGFR")
df <- df[df$sumstats %in% sumstats,]
df[df$sumstats=='chronic_kidney_disease',]$sumstats='CKD'
df[df$sumstats=='Microalbuminuria',]$sumstats='MAU'
df[df$sumstats=='MW_eGFR',]$sumstats='eGFR'

df$padj <- df$Coefficient_P_value
#df$padj <- p.adjust(df$Coefficient_P_value, method = "BH")
# plot each gwas trait result
for(gwas_trait in unique(df$sumstats)) {
  toplot <- dplyr::filter(df, sumstats == gwas_trait)

  # reorder the idents
  idents <- celltype
  toplot$celltype <- factor(toplot$celltype, levels = rev(idents))

  pval_threshold = -log10(0.05)

  # adjust pval and calculate pval threshold
  p1 <- ggplot(toplot, aes(x=celltype, y=-log10(padj), fill=ifelse(padj < 0.05,1,0))) + 
    geom_bar(stat="identity") + 
    geom_hline(yintercept = pval_threshold) +
    RotatedAxis() +
    ggtitle(gwas_trait) +
    theme_bw() +
    guides(fill = "none") +
    xlab("")+ylab("-log10(p_value)") +
    ylim(0,3.5) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text = element_text(size = 12),
          legend.title = element_blank()) +
    coord_flip()

  pdf(paste0(gwas_trait,".markerpeak.pdf"),width=4,height=5)
  print(p1)
  dev.off()
}



