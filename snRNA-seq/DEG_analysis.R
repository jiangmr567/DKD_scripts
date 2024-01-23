#########################################################################
# File Name: DEG_analysis.R
# 差异基因分析
#########################################################################
# load library
library(Seurat)
library(tibble)
library(dplyr)
library(future)
options(future.globals.maxSize = 30 * 1024^3)
plan(multisession, workers=16)

args = commandArgs(T)
outdir<-args[1]
setwd(outdir)

kidney <- readRDS("kidney_harmonyfinal.rds")
# 寻找DM与CTRL之间差异基因
Idents(kidney) <- "celltype"
lst <- list()
k=1
try(for(i in levels(kidney)){
	cellids <- subset(kidney,celltype==i)
	print(i)
	Idents(object = cellids) <- "state"
	lst[[k]] <- FindMarkers(cellids, ident.1 = "DM", ident.2 = "CTRL",logfc.threshold = 0)
	lst[[k]] <- lst[[k]] %>% rownames_to_column("gene")
	lst[[k]]$celltype <- i
	k=k+1
})
all_df=do.call(rbind,lst)
all_df = all_df %>% subset(p_val_adj<0.05 & abs(all_df$avg_log2FC)>0.25)
write.csv(all_df,"kidney_celltype_DM_CTRL_markers.csv")

# 寻找细胞类型之间的差异marker
data <- FindAllMarkers(kidney)
data = data %>% subset(p_val_adj<0.05 & data$avg_log2FC>0.25)
write.csv(data,"kidney_celltype_markers.csv")

deg_df <- read.csv("kidney_celltype_DM_CTRL_markers.csv",row.names=1)
head(deg_df)

###差异基因upset图
library(ComplexHeatmap)
deg_lst <- list()
for(cellname in unique(deg$celltype)){
	deg_lst[[cellname]] <- deg[deg$celltype==cellname,]$gene
}
mat <- list_to_matrix(deg_lst)

# cell-specific degs shared across cell types
m1 <- make_comb_mat(mat)
# filter by intersection size and set size
# only include intersections with a size of greater than 20, but include individual cell types
# identify combinations with only a single ident by summing all the occurrences of 1 in the string
df1 <- strsplit(names(comb_size(m1)), split="") %>%
  as.data.frame() 
df1 <- sapply(df1, as.numeric)
colsums <- colSums(df1)

# keep intersection set size of at least 15 and retain all groups with at least 1 deg
m1 <- m1[comb_size(m1) >= 15 | colsums == 1]

levels <- rev(c("Immune", "PODO", "ICB", "ICA", "Myofib", "Endo", "DCT", 
         "CNT", "PC", "PEC", "TAL", "PT_VCAM1","PT"))
ht <- draw(UpSet(m1, set_order=levels))
od <- column_order(ht)
cs <- comb_size(m1)
p3 <- draw(UpSet(m1, set_order=levels,top_annotation = upset_top_annotation(m1, add_numbers = TRUE)))

pdf("kidney_DEG_upset.pdf",height=6,width=8)
p3
dev.off()

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

###绘制差异基因火山图
library(ggplot2)
library(ggrepel)
celltype <- rev(levels)

#添加显著性标签
deg$celltype <- factor(deg$celltype,levels=celltype)
deg$label <- ifelse(deg$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")
#获取每个celltype中表达差异最显著的5个基因,去掉未知基因
deg1 <- deg[-grep("^ENSMFAG",deg$gene),]
deg1 %>% group_by(celltype) %>% top_n(n = 6, wt = abs(avg_log2FC)) -> top6_markers
#新增一列，将Top6的差异基因标记为2，其他的标记为1
deg$size <- case_when(!(deg$gene %in% top6_markers$gene)~ 1,
deg$gene %in% top6_markers$gene ~ 2)
#提取非Top6的基因表格
dt <- filter(deg,size==1)
#绘制每个Cluster Top6以外基因的散点火山图
p <- ggplot()+geom_jitter(data = dt,aes(x = celltype, y = avg_log2FC, color = label),size = 0.85,width =0.4)
#叠加每个Cluster Top6基因散点(将散点适当放大强调）
p <- ggplot()+scale_y_continuous(breaks=seq(-6, 7, 1))+
     geom_jitter(data = dt,aes(x = celltype, y = avg_log2FC, color = label),size = 0.85,width =0.4)+
     geom_jitter(data = top6_markers,aes(x = celltype, y = avg_log2FC, color = label),size = 1,width =0.4)
#根据图p中log2FC区间确定背景柱长度
dfbar<-data.frame(x=factor(celltype,levels=celltype),
                  y=c(2,2.75,3.45,2.35,1.85,2.15,3.45,2.05,2.55,1.9,2.55,2.65,2.5))
dfbar1<-data.frame(x=factor(celltype,levels=celltype),
                   y=c(-3.25,-2.1,-2.05,-2.6,-2.6,-2.9,-2.2,-2.15,-2.6,-2.1,-2.45,-3.4,-3.05))
#绘制背景柱
p1 <- ggplot()+geom_col(data = dfbar,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)+
      geom_col(data = dfbar1,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)
#把散点火山图叠加到背景柱上
p2 <- ggplot()+scale_y_continuous(breaks=seq(-6, 7, 1))+
      geom_col(data = dfbar,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)+
      geom_col(data = dfbar1,mapping = aes(x = x,y = y),fill = "#dcdcdc",alpha = 0.6)+
	  geom_jitter(data = dt,aes(x = celltype, y = avg_log2FC, color = label),size = 0.85,width =0.4)+
      geom_jitter(data = top6_markers,aes(x = celltype, y = avg_log2FC, color = label),size = 1,width =0.4)
#添加X轴的celltype色块标签：
dfcol<-data.frame(x=factor(celltype,levels=celltype),
                  y=0,label=celltype)
p3 <- p2 + geom_tile(data = dfcol,aes(x=x,y=y),height=0.4,color = "black",fill = col,alpha = 0.6,show.legend = F)
#给每个Cluster差异表达前Top6基因加上标签
p4 <- p3+geom_text_repel(data=top6_markers,aes(x=celltype,y=avg_log2FC,label=gene),force = 1.2,
      arrow = arrow(length = unit(0.008, "npc"),type = "open", ends = "last"))+
      scale_color_manual(name=NULL,values = c("red","black"))+labs(x="Celltype",y="avg_log2FC")+
      geom_text(data=dfcol,aes(x=x,y=y,label=label),size =5)
p5 <- p4+theme_minimal()+theme(axis.title = element_text(size = 14,color = "black",face = "bold"),
      axis.line.y = element_line(color = "black",size = 1.2),axis.line.x = element_blank(),axis.text.y=element_text(size=12),
	  axis.text.x = element_blank(),panel.grid = element_blank(),legend.position = "top",
	  legend.direction = "vertical",legend.justification = c(1,0),legend.text = element_text(size = 16))+
	  ggtitle("DM vs Control")+theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5))
ggsave("volcano_map.pdf",plot=p5,width=15,height=10)

deg$celltype <- factor(deg$celltype,levels=celltype)
deg_up <- subset(deg,avg_log2FC>0)
deg_down <- subset(deg,avg_log2FC<0)
table(deg_up$celltype)
table(deg_down$celltype)
###差异基因功能富集
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library(ggplot2)
library(stringr)
##功能注释
## 将SYMBOL转成ENTREZID
ids_up=bitr(deg_up$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
ids_down=bitr(deg_down$gene,'SYMBOL','ENTREZID','org.Hs.eg.db')
deg_up=merge(deg_up,ids_up,by.x='gene',by.y='SYMBOL')
gcSample=split(deg_up$ENTREZID, deg_up$celltype)
go_up<- compareCluster(gcSample,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05
)
deg_down=merge(deg_down,ids_down,by.x='gene',by.y='SYMBOL')
gcSample=split(deg_down$ENTREZID, deg_down$celltype)
go_down<- compareCluster(gcSample,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05
)
go_up  = setReadable(go_up, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
go_data_up<-as.data.frame(go_up[go_up@compareClusterResult$Count>=10,])
head(go_data_up)
#go_data_up <- subset(go_data_up,pvalue<0.01)
dim(go_data_up)
table(go_data_up$Cluster)
write.csv(go_data_up,"celltype_DEG_go_up_data.csv")
#画图
library(stringr)
p <- dotplot(go_up)+theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust = 0.5))+ggtitle("DEG UP Pathway Enrichment")+theme(plot.title = element_text(size = 16, face = "bold",hjust = 0.5))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))
ggsave("celltype_GO_up_plot.pdf",plot=p,width=12,height=12)

go_down  = setReadable(go_down, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
go_data_down<-as.data.frame(go_down)
go_data_down <- subset(go_data_down,p.adjust<0.05)
dim(go_data_down)
table(go_data_down$Cluster)
write.csv(go_data_down,"celltype_DEG_go_down_data.csv")
#画图
library(stringr)
p <- dotplot(go_down)+theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust = 0.5))+ggtitle("DEG DOWN Pathway Enrichment")+theme(plot.title = element_text(size = 16, face = "bold",hjust = 0.5))+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))
p
ggsave("celltype_GO_down_plot.pdf",plot=p,width=12,height=12)