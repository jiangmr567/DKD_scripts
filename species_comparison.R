#########################################################################
# File Name: species_comparison.R
# 物种比较
#########################################################################
# load library
library(ggpubr)
library(Matrix)
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ComplexHeatmap)
library(psych)
library(corrplot)
library(tidyverse)

args = commandArgs(T)
outdir<-args[1]
setwd(outdir)

kid <- readRDS("./data/human_mouse_macaca_kidney.rds")
head(kid)
table(kid$celltype1)
levels(kid)
## 人猴鼠RNA数据整合UMAP
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
colors <- colors[1:13]

kid$celltype1 <- factor(kid$celltype1,levels=rev(levels))
kid$state <- factor(kid$state,levels=c("human-CTRL", "monkey-CTRL", "mouse-CTRL", "human-DM","monkey-DM","mouse-DM"))
pdf('human_mouse_macaca_DimPlot.pdf',width = 12,height = 8)
DimPlot(kid, reduction = "umap", label = T,label.size =3,group.by = 'celltype1',cols=colors,pt.size=0.3)
DimPlot(kid, reduction = "umap", label = T,label.size =3, group.by = 'celltype1',pt.size=0.3,cols=colors,split.by = "state",ncol=3)
dev.off()

##三物种细胞类型相关性
#细胞类型markergene(p_val_adj<0.05 & avg_log2FC>0.25)
path <- "data/"
markers1 <- read.csv(paste0(path,"mouse_markers_Unite.csv"))
markers2 <- read.csv(paste0(path,"monkey_markers_Unite.csv"))
markers3 <- read.csv(paste0(path,"human_markers_Unite.csv"))
motop_all <- unique(markers1[markers1$avg_log2FC>0.25,]$gene)
matop_all <- unique(markers2[markers2$avg_log2FC>0.25,]$gene)
hutop_all <- unique(markers3[markers3$avg_log2FC>0.25,]$gene)
genes <- Reduce(intersect,list(motop_all,matop_all,hutop_all))   ###genes 2512

DefaultAssay(object = kid) <- "RNA"
Idents(kid) <- "cellanno1"
levels(kid)
levels(kid) <- c("PT-human","PT-monkey","PT-mouse","PT_VCAM1-human","PT_VCAM1-monkey","PT_VCAM1-mouse",
                 "TAL-human","TAL-monkey","TAL-mouse","PEC-human","PEC-monkey","PEC-mouse",
                 "ICA-human","ICA-monkey","ICA-mouse","ICB-human","ICB-monkey","ICB-mouse",
                 "PC-human","PC-monkey","PC-mouse","CNT-human","CNT-monkey","CNT-mouse",
                 "DCT-human","DCT-monkey","DCT-mouse","PODO-human","PODO-monkey","PODO-mouse",
                 "Endo-human","Endo-monkey","Endo-mouse","Myofib-human","Myofib-monkey","Myofib-mouse","Immune-human","Immune-monkey","Immune-mouse")

av1 <-AverageExpression(kid,group.by = "cellanno1",assays = "RNA")
av1=av1[[1]]


data_corr1 <- corr.test(av1[genes,levels(kid)], method="spearman", adjust="none")
data_r1 <- data_corr1$r          # 相关系数
data_p1 <- data_corr1$p          # p值
pmt <- data_r1
if (!is.null(pmt)){
  ssmt <- pmt>=0.9
  pmt[ssmt] <-'***'
  smt <- pmt >=0.75 & pmt <0.9
  pmt[smt] <- '**'
  mt <- pmt >=0.6 & pmt <0.75
  pmt[mt]<- '*'
  pmt[!ssmt&!smt&!mt]<- ''
} else {
  pmt <- F
}


celltype <- c()
species <- c()
for(i in levels(kid)){
  celltype <- c(celltype,strsplit(i,"-")[[1]][1])
  species <- c(species,strsplit(i,"-")[[1]][2])
}

annotation_col1 = data.frame(species=species,celltype = factor(celltype,levels=unique(celltype)))
rownames(annotation_col1) = levels(kid)
cell_col<-colors[1:length(levels(as.factor(annotation_col1$celltype)))]
names(cell_col)<-factor(unique(celltype),levels=unique(celltype))
ann_colors1 = list(celltype=cell_col,species = c(human = "#377EB8", monkey = "#6A3D9A",mouse="#33A02C"))

pdf('kidney_cor_analysis.pdf',width = 12,height = 10)
pheatmap(data_r1,name='corr',cluster_row = F,cluster_col = F, annotation_col = annotation_col1,annotation_colors = ann_colors1,
         fontsize = 12,display_numbers = pmt,main = "human vs monkey vs mouse correlation analysis")
dev.off()


library(UpSetR)
markers_h <- read.csv(paste0(path,"human_celltype_DM_CTRL_markers.csv"))
markers_ma <- read.csv(paste0(path,"monkey_celltype_DM_CTRL_markers.csv"))
markers_mo <- read.csv(paste0(path,"mouse_celltype_DM_CTRL_markers.csv"))
markers_h <- subset(markers_h,pct.1>0.2 | pct.2>0.2)
markers_ma <- subset(markers_ma,pct.1>0.2 | pct.2>0.2)
markers_mo <- subset(markers_mo,pct.1>0.2 | pct.2>0.2)

###疾病差异基因交并比
Idents(kid) <- "celltype1"
lst=list()
k=1
for(i in levels(kid)){
  a <- subset(markers_h,celltype==i)
  b <- subset(markers_ma,celltype==i)
  c <- subset(markers_mo,celltype==i)
  lst[[k]] <- data.frame(type=i)
  lst[[k]]$human_monkey <- length(intersect(a$gene,b$gene))/length(union(a$gene,b$gene))
  lst[[k]]$human_mouse <- length(intersect(a$gene,c$gene))/length(union(a$gene,c$gene))
  lst[[k]]$monkey_mouse <- length(intersect(b$gene,c$gene))/length(union(b$gene,c$gene))
  k=k+1
}
df <- as.data.frame(t(do.call(rbind,lst)),stringsAsFactors = F)[-1,]
colnames(df) <- levels(kid)
mtx <- apply(as.matrix(df),2,as.numeric)
rownames(mtx) <- c("human_monkey","human_mouse","monkey_mouse")

pdf('kidney_Jaccard_analysis.pdf',width = 8,height = 4)
corrplot(mtx,
         method="pie",col = COL1('YlOrRd'),
         tl.col="black", tl.srt=90,col.lim = c(0.3, 0.6), 
         is.corr = FALSE, addgrid.col = 'black' ,tl.cex = 1  
)
dev.off()

##PT细胞和内皮细胞GO分析
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library(ggplot2)
library(VennDiagram)
##PT共上调
a <- markers_h[(markers_h$avg_log2FC>0 & markers_h$celltype=='PT'),]
b <- markers_ma[(markers_ma$avg_log2FC>0 & markers_ma$celltype=='PT'),]
c <- markers_mo[(markers_mo$avg_log2FC>0 & markers_mo$celltype=='PT'),]
#组间交集元素获得
inter <- get.venn.partitions(list(a$gene,b$gene,c$gene))
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
inter[-c(5, 6)]$values[1]        ###哺乳动物共有基因
inter[-c(5, 6)]$values[5]        ###灵长类人猴特有基因
inter[-c(5, 6)]$values[4]        ###鼠特有基因
com_genes <- unlist(strsplit(inter[-c(5, 6)]$values[1],", "))
hu_ma_gene <- unlist(strsplit(inter[-c(5, 6)]$values[5],", "))
hu_mo_gene <- unlist(strsplit(inter[-c(5, 6)]$values[4],", "))
lst1=list()
lst1[[1]] <- data.frame(gene=com_genes,group="Mammals_PT_up")  
lst1[[2]] <- data.frame(gene=hu_ma_gene,group="Primates_PT_up")
lst1[[3]] <- data.frame(gene=hu_mo_gene,group="Rodents_PT_up")
deg_df <- do.call(rbind,lst1)
##PT共下调
a <- markers_h[(markers_h$avg_log2FC< 0 & markers_h$celltype=='PT'),]
b <- markers_ma[(markers_ma$avg_log2FC< 0 & markers_ma$celltype=='PT'),]
c <- markers_mo[(markers_mo$avg_log2FC< 0 & markers_mo$celltype=='PT'),]
#组间交集元素获得
inter <- get.venn.partitions(list(a$gene,b$gene,c$gene))
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
inter[-c(5, 6)]$values[1]        ###哺乳动物共有基因
inter[-c(5, 6)]$values[5]        ###灵长类人猴特有基因
inter[-c(5, 6)]$values[4]        ###鼠特有基因
com_genes <- unlist(strsplit(inter[-c(5, 6)]$values[1],", "))
hu_ma_gene <- unlist(strsplit(inter[-c(5, 6)]$values[5],", "))
hu_mo_gene <- unlist(strsplit(inter[-c(5, 6)]$values[4],", "))
lst2=list()
lst2[[1]] <- data.frame(gene=com_genes,group="Mammals_PT_down")  
lst2[[2]] <- data.frame(gene=hu_ma_gene,group="Primates_PT_down")
lst2[[3]] <- data.frame(gene=hu_mo_gene,group="Rodents_PT_down")
deg_df1 <- do.call(rbind,lst2)

##Endo共上调
a <- markers_h[(markers_h$avg_log2FC>0 & markers_h$celltype=='Endo'),]
b <- markers_ma[(markers_ma$avg_log2FC>0 & markers_ma$celltype=='Endo'),]
c <- markers_mo[(markers_mo$avg_log2FC>0 & markers_mo$celltype=='Endo'),]
#组间交集元素获得
inter <- get.venn.partitions(list(a$gene,b$gene,c$gene))
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
inter[-c(5, 6)]$values[1]        ###哺乳动物共有基因
inter[-c(5, 6)]$values[5]        ###灵长类人猴特有基因
inter[-c(5, 6)]$values[4]        ###鼠特有基因
com_genes <- unlist(strsplit(inter[-c(5, 6)]$values[1],", "))
hu_ma_gene <- unlist(strsplit(inter[-c(5, 6)]$values[5],", "))
hu_mo_gene <- unlist(strsplit(inter[-c(5, 6)]$values[4],", "))
lst1=list()
lst1[[1]] <- data.frame(gene=com_genes,group="Mammals_Endo_up")  
lst1[[2]] <- data.frame(gene=hu_ma_gene,group="Primates_Endo_up")
lst1[[3]] <- data.frame(gene=hu_mo_gene,group="Rodents_Endo_up")
deg_df2 <- do.call(rbind,lst1)
##Endo共下调
a <- markers_h[(markers_h$avg_log2FC< 0 & markers_h$celltype=='Endo'),]
b <- markers_ma[(markers_ma$avg_log2FC< 0 & markers_ma$celltype=='Endo'),]
c <- markers_mo[(markers_mo$avg_log2FC< 0 & markers_mo$celltype=='Endo'),]
#组间交集元素获得
inter <- get.venn.partitions(list(a$gene,b$gene,c$gene))
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
inter[-c(5, 6)]$values[1]        ###哺乳动物共有基因
inter[-c(5, 6)]$values[5]        ###灵长类人猴特有基因
inter[-c(5, 6)]$values[4]        ###鼠特有基因
com_genes <- unlist(strsplit(inter[-c(5, 6)]$values[1],", "))
hu_ma_gene <- unlist(strsplit(inter[-c(5, 6)]$values[5],", "))
hu_mo_gene <- unlist(strsplit(inter[-c(5, 6)]$values[4],", "))
lst2=list()
lst2[[1]] <- data.frame(gene=com_genes,group="Mammals_Endo_down")  
lst2[[2]] <- data.frame(gene=hu_ma_gene,group="Primates_Endo_down")
lst2[[3]] <- data.frame(gene=hu_mo_gene,group="Rodents_Endo_down")
deg_df3 <- do.call(rbind,lst2)

deg_df4 <- do.call(rbind,list(deg_df,deg_df1,deg_df2,deg_df3))
##功能注释
ids=bitr(deg_df4$gene,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID
deg_df4=merge(deg_df4,ids,by.x='gene',by.y='SYMBOL')
levels <- c("Mammals_Endo_down","Mammals_Endo_up","Primates_Endo_down","Primates_Endo_up",
            "Rodents_Endo_down","Rodents_Endo_up", "Mammals_PT_down", "Mammals_PT_up",
            "Primates_PT_down", "Primates_PT_up","Rodents_PT_down" ,  "Rodents_PT_up")
deg_df4$group <- factor(deg_df4$group,levels=levels)
##画上调的基因并进行功能注释
gcSample=split(deg_df4$ENTREZID, deg_df4$group)
go<- compareCluster(gcSample,
                    fun = "enrichGO",
                    OrgDb = "org.Hs.eg.db",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05
                    #qvalueCutoff = 0.05
)
go_data<-as.data.frame(go)
write.csv(go_data,"Endo_PT_go_data.csv")
#画图
p <- dotplot(go)+theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust = 0.5))+xlab("group")+
  scale_y_discrete(labels=function(x) str_wrap(x, width=60))
p
ggsave("Endo_PT_GO_plot.pdf",plot=p,width=12,height=12)

##
kid@meta.data$State <-"1"
kid@meta.data[grep("CTRL$",kid@meta.data$state),]$State <- "CTRL"
kid@meta.data[grep("DM$",kid@meta.data$state),]$State <- "DM"
kid@meta.data$group <- paste0(kid@meta.data$celltype1,'_',kid@meta.data$State)
kid@meta.data$celltype2 <- paste0(kid@meta.data$cellanno1,'_',kid@meta.data$State)

#基因表达小提琴图
PT <- subset(kid,celltype1=="PT")
gene3 <- c("PCK1","ALDOB","SLC16A9","INSR")
PT_count <- PT@assays$RNA@data[gene3,]
PT_count = as.data.frame(t(as.matrix(PT_count)))
PT_count$cellname = rownames(PT_count) 
head(PT_count)

PT_meta = PT@meta.data[,c("State","cellanno1")]
PT_meta$cellname = rownames(PT_meta)
PT_meta =merge(PT_count,PT_meta,by="cellname")
head(PT_meta)
PT_meta$celltype <- paste0(PT_meta$cellanno1,'-',PT_meta$State)

p <- ggviolin(PT_meta, "cellanno1", "SLC16A9", fill = "State",add = "boxplot",add.params = list(color = "black",width=0.08),
         ylim=c(0.25,6),width = 1.6) +
  theme_pubr(x.text.angle = 0, legend = "right") + 
  stat_pwc(aes(group = State), label = "p.signif",y.position =5.5) +
  xlab(NULL) + ggsci::scale_color_aaas() +
  ylab('Expression')+ggtitle("SLC16A9")

p1 <- ggviolin(PT_meta, "cellanno1", "PCK1", fill = "State",add = "boxplot",add.params = list(color = "black",width=0.08),
              ylim=c(0.25,6.5),width = 1.6) +
  theme_pubr(x.text.angle = 0, legend = "right") + 
  stat_pwc(aes(group = State), label = "p.signif",y.position =6) +
  xlab(NULL) + ggsci::scale_color_aaas() +
  ylab('Expression')+ggtitle("PCK1")

p2 <- ggviolin(PT_meta, "cellanno1", "ALDOB", fill = "State",add = "boxplot",add.params = list(color = "black",width=0.08),
              ylim=c(0.25,6),width = 1.6) +
  theme_pubr(x.text.angle = 0, legend = "right") + 
  stat_pwc(aes(group = State), label = "p.signif",y.position =5.5) +
  xlab(NULL) + ggsci::scale_color_aaas() +
  ylab('Expression')+ggtitle("ALDOB")

p3 <- ggviolin(PT_meta, "cellanno1", "INSR", fill = "State",add = "boxplot",add.params = list(color = "black",width=0.08),
               ylim=c(0.25,6),width = 1.6) +
  theme_pubr(x.text.angle = 0, legend = "right") + 
  stat_pwc(aes(group = State), label = "p.signif",y.position =5.5) +
  xlab(NULL) + ggsci::scale_color_aaas() +
  ylab('Expression')+ggtitle("INSR")

pdf("PT_maingene_ggviolin.pdf",width = 5,height = 4)
p
p1
p2
p3
dev.off()

###内皮细胞
Endo <- subset(kid,celltype1=="Endo")
gene3 <- c("CD74","TEK","DOCK4")
Endo_count <- Endo@assays$RNA@data[gene3,]
Endo_count = as.data.frame(t(as.matrix(Endo_count)))
Endo_count$cellname = rownames(Endo_count) 
head(Endo_count)

Endo_meta = Endo@meta.data[,c("State","cellanno1")]
Endo_meta$cellname = rownames(Endo_meta)
Endo_meta =merge(Endo_count,Endo_meta,by="cellname")
head(Endo_meta)
Endo_meta$celltype <- paste0(Endo_meta$cellanno1,'-',Endo_meta$State)

p <- ggviolin(Endo_meta, "cellanno1", "CD74", fill = "State",add = "boxplot",add.params = list(color = "black",width=0.08),
              ylim=c(0.25,6),width = 1.6) +
  theme_pubr(x.text.angle = 0, legend = "right") + 
  stat_pwc(aes(group = State), label = "p.signif",y.position =5.5) +
  xlab(NULL) + ggsci::scale_color_aaas() +
  ylab('Expression')+ggtitle("CD74")
p

p1 <- ggviolin(Endo_meta, "cellanno1", "DOCK4", fill = "State",add = "boxplot",add.params = list(color = "black",width=0.08),
               ylim=c(0.25,6.5),width = 1.6) +
  theme_pubr(x.text.angle = 0, legend = "right") + 
  stat_pwc(aes(group = State), label = "p.signif",y.position =6) +
  xlab(NULL) + ggsci::scale_color_aaas() +
  ylab('Expression')+ggtitle("DOCK4")
p1

p2 <- ggviolin(Endo_meta, "cellanno1", "TEK", fill = "State",add = "boxplot",add.params = list(color = "black",width=0.08),
               ylim=c(0.25,6),width = 1.6) +
  theme_pubr(x.text.angle = 0, legend = "right") + 
  stat_pwc(aes(group = State), label = "p.signif",y.position =5.5) +
  xlab(NULL) + ggsci::scale_color_aaas() +
  ylab('Expression')+ggtitle("TEK")
p2

pdf("Endo_maingene_ggviolin.pdf",width = 5,height = 4)
p
p1
p2
dev.off()

###重置点图--共有基因INSR
head(DotPlot(kid, features = 'INSR',group.by = 'celltype2')$data)
df <- DotPlot(kid, features = 'INSR',group.by = 'celltype2')$data
celltype <- c()
state <- c()
species <- c()
for(i in df$id){
  celltype <- c(celltype,strsplit(i,"-")[[1]][1])
  state <- c(state,strsplit(i,"-")[[1]][2])
}
df$celltype <- celltype
level <- c("PT","PT_VCAM1","TAL","ICA","ICB","PEC","PC","CNT","DCT","PODO","Endo","Myofib","Immune")
df$celltype <- factor(df$celltype,levels = rev(level))
df$state <- state
for(i in df$state){
  species <- c(species,strsplit(i,"_")[[1]][1])
}
df$species <- species
p <- ggplot(df,aes(x=state,y =celltype,size = pct.exp, color = avg.exp.scaled))+
  geom_point() + 
  scale_size("Percent\nExpressed", range = c(0,10)) + #调整绘图点的相对大小
  scale_color_gradientn(colours = viridis::viridis(30),
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "Average\nExpression") +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("") + theme_bw() +
  #scale_y_continuous(breaks = 1:length(levels(df$celltype)),labels = levels(df$celltype),sec.axis = dup_axis())+ #复制 y轴 代替边框效果
  facet_grid(~species, scales="free_x",space = "free")+theme_classic() +
  theme(
    axis.text.x = element_text(size=12, color="black",face="bold"),#x轴标签样式
    axis.text.y = element_text(size=12, color="darkblue",face="bold"),
    axis.ticks.y = element_blank(),#坐标轴刻度
    axis.text.y.right = element_blank(),#坐标轴标签隐藏
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = 'grey30',size = 0.2), #坐标轴轴线样式
    
    panel.spacing=unit(0, "mm"), #分面图图块间距
    strip.text.x = element_text(size=15, face="bold",color = "#FFFFFF",
                                vjust = 0.5,margin = margin(b = 3,t=3)),#分面标签样式
    strip.background = element_rect(colour="grey30", fill="grey60",size = 1)
  )+RotatedAxis()
g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
cols <- colors[1:3]

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- cols[i]
}
pdf("common_gene_INSR_dot.pdf",width=5,height=6)
plot(g)
dev.off()

###重置点图--共有基因NR3C1
head(DotPlot(kid, features = 'NR3C1',group.by = 'celltype2')$data)
df <- DotPlot(kid, features = 'NR3C1',group.by = 'celltype2')$data
celltype <- c()
state <- c()
species <- c()
for(i in df$id){
  celltype <- c(celltype,strsplit(i,"-")[[1]][1])
  state <- c(state,strsplit(i,"-")[[1]][2])
}
df$celltype <- celltype
level <- c("PT","PT_VCAM1","TAL","ICA","ICB","PEC","PC","CNT","DCT","PODO","Endo","Myofib","Immune")
df$celltype <- factor(df$celltype,levels = rev(level))
df$state <- state
for(i in df$state){
  species <- c(species,strsplit(i,"_")[[1]][1])
}
df$species <- species
p <- ggplot(df,aes(x=state,y =celltype,size = pct.exp, color = avg.exp.scaled))+
  geom_point() + 
  scale_size("Percent\nExpressed", range = c(0,10)) + #调整绘图点的相对大小
  scale_color_gradientn(colours = viridis::viridis(30),
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "Average\nExpression") +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("") + theme_bw() +
  #scale_y_continuous(breaks = 1:length(levels(df$celltype)),labels = levels(df$celltype),sec.axis = dup_axis())+ #复制 y轴 代替边框效果
  facet_grid(~species, scales="free_x",space = "free")+theme_classic() +
  theme(
    axis.text.x = element_text(size=12, color="black",face="bold"),#x轴标签样式
    axis.text.y = element_text(size=12, color="darkblue",face="bold"),
    axis.ticks.y = element_blank(),#坐标轴刻度
    axis.text.y.right = element_blank(),#坐标轴标签隐藏
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = 'grey30',size = 0.2), #坐标轴轴线样式
    
    panel.spacing=unit(0, "mm"), #分面图图块间距
    strip.text.x = element_text(size=15, face="bold",color = "#FFFFFF",
                                vjust = 0.5,margin = margin(b = 3,t=3)),#分面标签样式
    strip.background = element_rect(colour="grey30", fill="grey60",size = 1)
  )+RotatedAxis()
g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
cols <- colors[1:3]

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- cols[i]
}
pdf("common_gene_NR3C1_dot.pdf",width=5,height=6)
plot(g)
dev.off()