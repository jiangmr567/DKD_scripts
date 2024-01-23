#########################################################################
# File Name: monkey2human.R
# 将猴的基因组坐标转换为hg19基因组坐标
#########################################################################
# load library
library(Seurat)
library(dplyr)
library(data.table)
library(GenomicFeatures)
library(parallel)
library(clusterProfiler)
#library(tidyverse)
library(patchwork)

###猴与人类同源gene转换
library(biomaRt)
library(dplyr)
library(stringr)

args = commandArgs(T)
outdir<-args[1]
rna_path <- args[2]
atac_path <- args[3]
setwd(outdir)

### 获取数据
## Nov 2020 v102 基因组版本为Macaca_fascicularis_5.0
human = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "nov2020.archive.ensembl.org")
cyno = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mfascicularis_gene_ensembl",  host = "nov2020.archive.ensembl.org")

## 构建转化函数
human2cyno <- function(x){
  genes_cyno = getLDS(attributes = c("hgnc_symbol"), 
                      filters = "hgnc_symbol", 
                      values = x , mart = human, 
                      attributesL = c("external_gene_name"), 
                      martL = cyno, uniqueRows=T
  )
  colnames(genes_cyno) = c("human_gene", "gene")
  return(genes_cyno)
}

###使用markergene卡值选择：avg_log2FC>1 & p_val_adj<0.05
DEG_df <- read.csv(paste0(rna_path,"/kidney_celltype_markers.csv"),row.names = 1)
genes <- DEG_df$gene
gene_match <- human2cyno(genes)
head(gene_match)
head(DEG_df)
DEG_df <- na.omit(merge(DEG_df,gene_match,by='gene'))
DEG_df1 <- subset(DEG_df,p_val_adj<0.05 & avg_log2FC>1)
write.csv(DEG_df1,"./00.gene2human/celltype_marker2human.csv")

#按细胞类型输出每种细胞类型的marker,将每种细胞类型的marker转成chr.start.end三列
DEG_df <- read.csv("./00.gene2human/celltype_DEG2human.csv",row.names=1)  ###已经转成人同源基因
tss <- read.table("./00.gene2human/hg19_tss.bed.gz")
deg_df=merge(tss, DEG_df, by.x='V3', by.y='human_gene')
deg_df$chr  = paste0('chr',deg_df$V1)
deg_df$start= deg_df$V2-2000
deg_df$end  = deg_df$V2+2000
deg_df=split(deg_df,deg_df$celltype)
lapply(names(deg_df),function(i)write.table(deg_df[[i]][,c('chr','start','end')],row.names=F,col.names=F,quote=F,sep="\t",file=paste0("01.celtype.DEG.star.end/celltype_DEG/",i,'.bed')))

##使用markerpeak去转成人的坐标(FDR<0.01 & Log2FC>4)
file <- paste0(atac_path,"/markerPeak_celltype_final/")
ff <- list.files(file,"_markerPeak.bed$")
peak_lst <- list()
for(i in seq_along(ff)){
	peak_lst[[i]] <- read.table(paste0(file,ff[i]),sep='\t')
	peak_lst[[i]]$celltype <- gsub('_markerPeak.bed',"",ff[i])
}
peak_df <- do.call(rbind,peak_lst)
colnames(peak_df) <- c("seqnames","start","end","cellname","Log2FC","FDR","celltype")
peak_df <- subset(peak_df,FDR<0.01 & Log2FC>4)
write.table(peak_df, file="./00.gene2human/celltype_markerpeak.bed",row.names=F,col.names=F,sep="\t",quote=F)

##利用https://genome.ucsc.edu/cgi-bin/hgLiftOver 将peak文件转换成人类h19坐标，读入转出结果
peak_df <- read.table("./00.gene2human/celltype_markerpeak_h19.bed",sep='\t')
peak_df <- peak_df[,c(1,2,3,7)]
colnames(peak_df) <- c("seqnames","start","end","celltype")
for (ct in unique(peak_df$celltype)){
	print(ct)
	df = peak_df %>% filter(celltype %in% ct)
	peak.data= df[,c("seqnames","start","end")]
	ct=gsub(" ","_",ct)
	ct=gsub("/","_",ct)
	ct=gsub("-","_",ct)
	write.table(peak.data,paste("01.celtype.DEG.star.end/celltype_peak/",ct,".bed",sep=""),quote=F,row.names=F,col.names=F,sep="\t")
}