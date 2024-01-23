#########################################################################
# File Name: cellchat.R
# 细胞相互作用分析
#########################################################################
# load library
library(CellChat)
library(svglite)
library(Seurat)
library(patchwork) 
#library(reshape)
library(ComplexHeatmap)
set.seed(123) 

args = commandArgs(T)
outdir<-args[1]
setwd(outdir)

rna <- readRDS("kidney_harmonyfinal.rds")
table(rna$state)
###将疾病对照细胞数控制到同一纬度
Idents(rna) <- "state"
rna <- subset(rna,downsample=10000)
table(rna$state)
kidney <- rna
Idents(kidney) <-'celltype'
kidney_DM <- subset(kidney,subset=state=="DM")
kidney_CTRL <- subset(kidney,subset=state=="CTRL")
state <- c(kidney_DM,kidney_CTRL)

k=1
object.list <- list()
for(i in state){
	data.input <- GetAssayData(i, assay = "RNA", slot = "data") # normalized data matrix
	labels <- Idents(i)
	meta <- data.frame(group = labels, row.names = names(labels)) 

	# 构建cellChat对象
	cellchat <- createCellChat(object = data.input)
	# 细胞类型信息加入object对象中
	cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
	cellchat <- setIdent(cellchat, ident.use = "labels")
	# 导入配受体数据库
	CellChatDB <- CellChatDB.human

	# 选取分泌信号通路进行下游分析
	#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
	#cellchat@DB <- CellChatDB.use    # set the used database in the object
	# 选取全部信号通路进行下游分析
	cellchat@DB <- CellChatDB
	unique(CellChatDB$interaction$annotation)
	# 首先在一个细胞群中识别过表达的配体或受体，然后识别其相互作用；
	# 可将基因表达数据投射到蛋白-蛋白相互作用(PPI)网络上。

	cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost 这一步是必须的
	#future::plan("multiprocess", workers = 8) # do parallel  
	# 识别过表达基因
	cellchat <- identifyOverExpressedGenes(cellchat)
	# 识别过表达基因的互作
	cellchat <- identifyOverExpressedInteractions(cellchat)

	# 通过PPI对基因表达进行smooth
	# cellchat <- projectData(cellchat, PPI.human)  # 注意如果进行PPI smooth，则进行computeCommunProb计算时，raw.use = FALSE

	# CellChat通过给每个相互作用赋一个概率值并进行排列、测试，来推断具有生物意义的细胞通讯。
	# 通过置换检验来推断生物意义上的细胞-细胞通信。

	# 计算通讯概率
	cellchat <- computeCommunProb(cellchat, raw.use = TRUE)  # 默认 raw.use = TRUE

	# 对细胞进行过滤
	cellchat <- filterCommunication(cellchat, min.cells = 5)

	# 通过计算的每个配体-受体通信概率来计算信号通路层面的通信概率
	# 每个信号通路分别存储在“net”和“netP”槽中
	cellchat <- computeCommunProbPathway(cellchat)

	# 计算聚集的细胞通讯网络
	cellchat <- aggregateNet(cellchat)
	cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
	object.list[k] <- cellchat
	k=k+1
}
saveRDS(object.list[[1]],"./cellchat/CTRL_cellchat.rds")
saveRDS(object.list[[2]],"./cellchat/DM_cellchat.rds")

##不同生物条件细胞间通讯网络的相互作用总数和相互作用强度
names(object.list) <- c('CTRL','DM')
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("./cellchat/state_cellchat_strength.pdf")
gg1+gg2
dev.off()

saveRDS(cellchat, file = "./cellchat/kidney_cellchat_comparisonAnalysis.rds")

###相互作用权重
weight_ctrl <- cellchat@net$CTRL$weight
weight_dm <- cellchat@net$DM$weight
weight_ctrl
weight_dm
p1 <- Heatmap(as.matrix(weight_ctrl),name="weight",
         col=c("#4264a3","white","#9b1a28"),
         cluster_rows = F,
         cluster_columns = F,
         show_column_names = T,
         show_row_names = T,
         row_title="Sources",
         column_title = "Interaction weights of CTRL",
         row_names_side = "left",column_names_side = "bottom",
         column_names_max_height = unit(200, "cm"))
         
pdf("./cellchat/cellchat_weight_CTRL_heatmap.pdf",width = 5.5,height = 5)
p1
dev.off()
p2 <- Heatmap(as.matrix(weight_dm),name="weight",
         col=c("#4264a3","white","#9b1a28"),
         cluster_rows = F,
         cluster_columns = F,
         show_column_names = T,
         show_row_names = T,
         row_title="Sources",
         column_title = "Interaction weights of DM",
         row_names_side = "left",column_names_side = "bottom",
         column_names_max_height = unit(200, "cm"))

pdf("./cellchat/cellchat_weight_DM_heatmap.pdf",width = 5.5,height = 5)
p2
dev.off()

##在不同状态之间发送或接收信号发生显着变化的细胞群
num.link <- unlist(sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)}))
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
pdf("./cellchat/celltype_outin_change.pdf",height=5,width=6)
wrap_plots(plots = gg)
dev.off()

##识别上调和下调的信号配体-受体对
pos.dataset = "DM"
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
net <- netMappingDEG(cellchat, features.name = features.name)
net.all <- subsetCommunication(cellchat, net = net)
net_filt <- subset(net.all,ligand %in% c("VEGFA","COL4A3","COL4A4","COL4A5") & receptor %in% c("FLT1_KDR","FLT1","KDR","ITGA2_ITGB1","ITGAV_ITGB8","ITGA3_ITGB1"))
pairLR.use.filt = net_filt[, "interaction_name", drop = F]
pdf("./cellchat/PODO_Endo_regulated_signaling.pdf",width=4,height=5)
netVisual_bubble(cellchat, pairLR.use = pairLR.use.filt, targets.use =c(8,12), sources.use = c(12), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = "all-regulated signaling in DM")
dev.off()

## 配受体基因表达差异
seu <- subset(rna,celltype %in% c("Endo","PODO"))
seu$group <- paste0(seu$celltype,'-',seu$state)
genes <- c("FLT1","KDR","VEGFA","ITGA2","ITGA3","ITGB1","ITGAV","ITGB8","COL4A5","COL4A4","COL4A3")
head(DotPlot(seu, features = genes,group.by = 'group')$data)
df <- DotPlot(seu, features = genes,group.by = 'group')$data
df$gene <- rownames(df)
celltype <- c()
state <- c()
for(i in df$id){
  celltype <- c(celltype,strsplit(i,"-")[[1]][1])
  state <- c(state,strsplit(i,"-")[[1]][2])
}
df$gene <- df$features.plot
df$state <- state
df$celltype <- celltype
df$celltype <- factor(df$celltype,levels=c("Endo","PODO"))
head(df)

## Sequential colormaps:
solarExtra = c('#3361A5', '#248AF3', '#14B3FF', '#88CEEF', '#C1D5DC', 
            '#EAD397', '#FDB31A', '#E42A2A', '#A31D1D')#buencolors
            
p <- ggplot(df,aes(x=state,y =gene,size = pct.exp, color = avg.exp.scaled))+
  geom_point() + 
  scale_size("Percent\nExpressed", range = c(0,6)) + #调整绘图点的相对大小
  scale_color_gradientn(colours = solarExtra,
                        guide = guide_colorbar(ticks.colour = "black",frame.colour = "black"),
                        name = "Average\nExpression") +
  cowplot::theme_cowplot() + 
  ylab("") + xlab("") + theme_bw() +
  facet_grid(~celltype, scales="free_x",space = "free") +
  theme(
    axis.text.x = element_text(size=10, color="black",face="bold"),#x轴标签样式
    axis.text.y = element_text(size=10, color="darkblue",face="bold"),
    axis.ticks.y = element_blank(),#坐标轴刻度
    axis.text.y.right = element_blank(),#坐标轴标签隐藏
    axis.ticks.x = element_blank(),
    axis.line = element_line(colour = 'grey30',size = 0.2), #坐标轴轴线样式
    
    panel.spacing=unit(0, "mm"), #分面图图块间距
    strip.text.x = element_text(size=12, face="bold",color = "#FFFFFF",
                                vjust = 0.5,margin = margin(b = 3,t=3)),#分面标签样式
    strip.background = element_rect(colour="grey30", fill="grey60",size = 1)
  )+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())
g <- ggplot_gtable(ggplot_build(p))
strips <- which(grepl('strip-', g$layout$name))
cols <- c("#BF5B17","#E6AB02")

for (i in seq_along(strips)) {
  k <- which(grepl('rect', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  l <- which(grepl('titleGrob', g$grobs[[strips[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strips[i]]]$grobs[[1]]$children[[k]]$gp$fill <- cols[i]
}

pdf("./cellchat/PODO_Endo_regulated_signaling_LR_dotplot.pdf",width=4,height=5)
plot(g)
dev.off()