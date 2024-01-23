#########################################################################
# File Name: SoupX.R
# 降低环境RNA的噪音
#########################################################################
# load library
library(SoupX)
library(Seurat)
set.seed(123) 

indir="AUTO_file"
outdir="soupx_result"

setwd(outdir)
## 参数简介
# toc是分析矩阵，即有过滤的矩阵
# tod是全矩阵，即没有任何过滤的矩阵
# rho是污染比例系数，可自行设置，如果不设置则会自动计算
for(file in list.files(indir)){
	file1 <- list.files(paste0(indir,"/",file))
    print(file1)
	toc <- Read10X(paste0(indir,"/",file,"/",file1,"/04.Matrix"),gene.column=1)
	print(length(unique(rownames(toc))))
	tod <- Read10X(paste0(outdir,"/",file,"/ALL"),gene.column=1)
	print(length(unique(rownames(tod))))
    tod <- tod[rownames(toc),]      # 取toc里面的所有基因，保证基因名一致
	# SoupX帮助文档建议提供分析矩阵的聚类亚群分组，因此这里利用分析矩阵做一个简单聚类
	all <- toc
	all <- CreateSeuratObject(all)
	all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
	all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 3000)
	all.genes <- rownames(all)
	all <- ScaleData(all, features = all.genes)

	all <- RunPCA(all, features = VariableFeatures(all), npcs = 40, verbose = F)
	all <- FindNeighbors(all, dims = 1:30)
	all <- FindClusters(all, resolution = 0.5)
	all <- RunUMAP(all, dims = 1:30)

	matx <- all@meta.data
	sc = SoupChannel(tod, toc)
	sc = setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))
	# 自动计算污染比例系数
	tryCatch(
	{sc = autoEstCont(sc)}, 
	error=function(e) {
		#因为自动计算经常会报错，所以如果报错则设置rho为0.2
		sc <<- setContaminationFraction(sc, 0.2)        # 设置为全局变量
		print("autoEstCont Error !")})
	out = adjustCounts(sc)
	saveRDS(sc,paste0(file,"/sc_",file,".rds"))
	# 保存校正后的矩阵，输出为10X格式
	DropletUtils:::write10xCounts(paste0(file,"/soupX_matrix"), out,version="3",overwrite = T)
}