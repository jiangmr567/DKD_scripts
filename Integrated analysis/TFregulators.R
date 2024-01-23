#########################################################################
# File Name: TFregulators.R
# TF调节因子潜在调控靶基因的鉴定
#########################################################################
# load library
library(ArchR)
library(parallel)
library(Seurat)
library(pryr) 
library(cowplot)
addArchRThreads(threads = 16)
set.seed(123)
library("BSgenome.Mfascicularis.NCBI.5.0")
assignInNamespace(".h5read", function (file = NULL, name = NULL, method = "fast", index = NULL, 
    start = NULL, block = NULL, count = NULL) {
    res <- rhdf5::h5read(file = file, name = name, index = index, 
            start = start, block = block, count = count)
    o <- h5closeAll()
    return(res)
}, ns="ArchR")
seqnames(BSgenome.Mfascicularis.NCBI.5.0)=gsub("MFA","chr",seqnames(BSgenome.Mfascicularis.NCBI.5.0))

args = commandArgs(T)
outdir<-args[1]
setwd(outdir)

genomeAnnotation <- readRDS("../data/Macaca_genomeAnnotation")
geneAnnotation <- readRDS("../data/Macaca_geneAnnotation")

proj <- readRDS("Save-ArchR-Project.rds")

###参考https://github.com/GreenleafLab/scScalpChromatin/blob/main/Figure_3_TFregulators.R
source("../utils/plotting_config.R")
source("../utils/misc_helpers.R")
source("../utils/matrix_helpers.R")
source("../utils/archr_helpers.R")
source("../utils/GO_wrappers.R")

pointSize <- 1.0
corrCutoff <- 0.45 # Used in labeling peak2gene links
dir.create("TFregulators")

# GeneIntegration Matrix: (rows gene names x cols cell names)
GIMatrix <- getMatrixFromProject(proj, useMatrix="GeneIntegrationMatrix")
GImat <- assays(GIMatrix)$GeneIntegrationMatrix
rownames(GImat) <- rowData(GIMatrix)$name
GImat <- as(GImat[Matrix::rowSums(GImat) > 0,], "sparseMatrix") # Remove unexpressed genes

regulators <- c('HNF4A_582','HNF4G_583','HNF1A_85','HNF1B_86','PPARD_436',"RXRG_265","ETV1_568","ETV4_569","ELF1_565","IKZF1_397","GABPA_576")
regulators

deviationsMatrix <- deviationsMatrix[regulators,]
# Identify pseudobulks for performing matrix correlations
knn_groups <- getLowOverlapAggregates(proj, target.agg=500, k=100, 
  overlapCutoff=0.8, dimReduc="Harmony")

kgrps <- unique(knn_groups$group)
# GeneIntegrationMatrix
GIMatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowMeans(GImat[,use_cells])
  }) %>% do.call(cbind,.)
colnames(GIMatPsB) <- kgrps

# In rare instances, we can get pseudo-bulked genes that have zero averages
GIMatPsB <- GIMatPsB[Matrix::rowSums(GIMatPsB) > 0,]

# DeviationsMatrix
DevMatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowMeans(deviationsMatrix[,use_cells])
  }) %>% do.call(cbind,.)
colnames(DevMatPsB) <- kgrps

# Perform chromVAR deviations to Integrated RNA correlation analysis:
start <- Sys.time()
geneCorMat <- cor2Matrices(DevMatPsB, GIMatPsB)
colnames(geneCorMat) <- c("motifName", "symbol", "Correlation", "FDR")
end <- Sys.time()
message(sprintf("Finished correlations in %s minutes.", round((end  - start)/60.0, 2)))

allGenes <- rownames(GIMatPsB) %>% sort() # Already filtered to only expressed genes

# Get locations of motifs of interest:
motifPositions <- getPositions(proj, name="Motif")
motifGR <- stack(motifPositions, index.var="motifName")

# Get peak to gene GR
p2gGR <- getP2G_GR(proj, corrCutoff=corrCutoff)

calculateLinkageScore <- function(motifLocs, p2gGR){
  # Calculate Linkage Score (LS) for each gene in p2gGR with regards to a motif location GR
  ###################################
  # For a given gene, the LS = sum(corr peak R2 * motifScore)
  ol <- findOverlaps(motifLocs, p2gGR, maxgap=0, type=c("any"), ignore.strand=TRUE)
  olGenes <- p2gGR[to(ol)]
  olGenes$motifScore <- motifLocs[from(ol)]$score
  olGenes$R2 <- olGenes$Correlation**2 # All p2g links here are already filtered to only be positively correlated
  LSdf <- mcols(olGenes) %>% as.data.frame() %>% group_by(symbol) %>% summarise(LS=sum(R2*motifScore)) %>% as.data.frame()
  LSdf <- LSdf[order(LSdf$LS, decreasing=TRUE),]
  LSdf$rank <- 1:nrow(LSdf)
  return(LSdf)
}

calculateMotifEnrichment <- function(motifLocs, p2gGR){
  # Calculate Motif enrichment per gene
  ###################################
  # For a given gene, calculate the hypergeometric enrichment of motifs in 
  # linked peaks (generally will be underpowered)
  motifP2G <- p2gGR[overlapsAny(p2gGR, motifLocs, maxgap=0, type=c("any"), ignore.strand=TRUE)]
  m <- length(motifP2G) # Number of possible successes in background
  n <- length(p2gGR) - m # Number of non-successes in background

  motifLinks <- motifP2G$symbol %>% getFreqs()
  allLinks <- p2gGR$symbol %>% getFreqs()
  df <- data.frame(allLinks, motifLinks=motifLinks[names(allLinks)])
  df$motifLinks[is.na(df$motifLinks)] <- 0
  df$mLog10pval <- apply(df, 1, function(x) -phyper(x[2]-1, m, n, x[1], lower.tail=FALSE, log.p=TRUE)/log(10))
  df <- df[order(df$mLog10pval, decreasing=TRUE),]
  df$symbol <- rownames(df)
  return(df)
}

##############近端小管PT################
library(msigdbr)
library(dplyr)
db <- msigdbr()
Glycolysis <- db %>% filter(gs_name == "KEGG_GLYCOLYSIS_GLUCONEOGENESIS") %>% pull(human_gene_symbol) 
Oxidative <- db %>% filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") %>% pull(human_gene_symbol)
Metabolism <- db %>% filter(gs_name == "HALLMARK_FATTY_ACID_METABOLISM") %>% pull(human_gene_symbol)

deg <- read.csv("../RNA/result/kidney_celltype_DM_CTRL_markers.csv",row.names=1)
deg <- subset(deg,abs(avg_log2FC)>0.5 & p_val_adj<0.05)
pt_deg <- subset(deg,celltype=='PT')
head(pt_deg)
markerGenes <- intersect(c(Glycolysis, Oxidative,Metabolism),pt_deg$gene) %>% unique() %>% sort()
markerGenes <- c(markerGenes,"SLC34A1","SLC13A3","LRP2","SLC17A3","SLC4A4") %>% unique()
markerGenes

regulators <- c('HNF4A_582','HNF4G_583','HNF1A_85','HNF1B_86','PPARD_436',"RXRG_265")
# Store results for each TF
res_list <- list()

##########################################################################################
for(motif in regulators){
  motif_short <- strsplit(motif,"_")[[1]][1]
  # First get motif positions
  motifLocs <- motifGR[motifGR$motifName == motif]
  # Calculate Linkage Score for motif
  LS <- calculateLinkageScore(motifLocs, p2gGR)
  # Get just genes correlated to motif
  motifGeneCorDF <- geneCorMat[geneCorMat$motifName == motif,]
  plot_df <- merge(LS, motifGeneCorDF, by="symbol", all.x=TRUE)
  # Calculate motif enrichment per gene
  ME <- calculateMotifEnrichment(motifLocs, p2gGR)
  plot_df <- merge(plot_df, ME, by="symbol", all.x=TRUE)
  plot_df <- plot_df[,c("symbol", "LS", "Correlation", "FDR", "mLog10pval")]
  plot_df$toLabel <- "NO"
  topN <- 5
  plot_df <- plot_df[order(plot_df$LS, decreasing=TRUE),]
  plot_df$rank_LS <- 1:nrow(plot_df)
  plot_df$toLabel[1:topN] <- "YES"
  plot_df <- plot_df[order(plot_df$Correlation, decreasing=TRUE),]
  plot_df$rank_Corr <- 1:nrow(plot_df)
  plot_df$toLabel[1:topN] <- "YES"
  plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=TRUE),]
  plot_df$rank_Pval <- 1:nrow(plot_df)
  plot_df$toLabel[1:10] <- "YES"
  plot_df$meanRank <- apply(plot_df[,c("rank_LS", "rank_Corr", "rank_Pval")], 1, mean)
  plot_df <- plot_df[order(plot_df$meanRank, decreasing=FALSE),]
  plot_df$toLabel[1:topN] <- "YES"
  # Label any marker genes in window of interest
  LS_window <- quantile(plot_df$LS, 0.8)
  corr_window <- 0.25
  pos_top_genes <- plot_df[plot_df$LS > LS_window & plot_df$Correlation > corr_window,]$symbol
  write.csv(pos_top_genes, paste0('./TFregulators', sprintf("/regulatory_targets_%s.csv", motif_short)))
  #neg_top_genes <- plot_df[plot_df$LS > LS_window & -plot_df$Correlation > corr_window,]$symbol
  if(nrow(plot_df[plot_df$symbol %in% pos_top_genes & plot_df$symbol %in% markerGenes,]) > 0){
    plot_df[plot_df$symbol %in% pos_top_genes & plot_df$symbol %in% markerGenes,]$toLabel <- "YES"
  }
  res_list[[motif_short]] <- pos_top_genes # Save regulatory targets
  # Save dataframe of results
  save_df <- plot_df[plot_df$symbol %in% pos_top_genes,c(1:5)]
  save_df <- save_df[order(save_df$Correlation, decreasing=TRUE),]
  saveRDS(save_df, paste0('./TFregulators', sprintf("/regulatory_targets_%s.rds", motif_short)))
  plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=FALSE),]
  # Label motif as well
  plot_df$toLabel[which(plot_df$symbol == motif_short)] <- "YES"
  plot_df$symbol[which(plot_df$toLabel == "NO")] <- ""
  # Threshold pvalue for plotting
  maxPval <- 5
  plot_df$mLog10pval <- ifelse(plot_df$mLog10pval > maxPval, maxPval, plot_df$mLog10pval)
  #Plot results
  p <- (
    ggplot(plot_df, aes(x=Correlation, y=LS, color=mLog10pval)) 
      #+ geom_point(size = 2)
      + ggrastr::geom_point_rast(size=2)
      + ggrepel::geom_text_repel(
          data=plot_df[plot_df$toLabel=="YES",], aes(x=Correlation, y=LS, label=symbol), 
          #data = plot_df, aes(x=Correlation, y=LS, label=symbol), #(don't do this, or the file will still be huge...)
          size=3,
          point.padding=0, # additional pading around each point
          box.padding=0.5,
          min.segment.length=0, # draw all line segments
          max.overlaps=Inf, # draw all labels
          #nudge_x = 2,
          color="black"
      ) 
      + geom_vline(xintercept=0, lty="dashed") 
      + geom_vline(xintercept=corr_window, lty="dashed", color="red")
      + geom_hline(yintercept=LS_window, lty="dashed", color="red")
      + theme_BOR(border=FALSE)
      + theme(panel.grid.major=element_blank(), 
              panel.grid.minor= element_blank(), 
              plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
              aspect.ratio=1.0,
              #legend.position = "none", # Remove legend
              axis.text.x = element_text(angle=90, hjust=1))
      + ylab("Linkage Score") 
      + xlab("Motif Correlation to Gene") 
      + scale_color_gradientn(colors=cmaps_BOR$zissou, limits=c(0, maxPval))
      + scale_y_continuous(expand = expansion(mult=c(0,0.05)))
      + scale_x_continuous(limits = c(-0.85, 0.955)) # Force plot limits
      + ggtitle(sprintf("%s putative targets", motif_short))
      )
  # Positively regulated genes:
  upGO <- calcTopGo(allGenes, interestingGenes=pos_top_genes, nodeSize=5, ontology="BP")
  upGO <- upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]
  up_go_plot <- topGObarPlot(upGO, cmap=cmaps_BOR$comet, nterms=3, border_color="black", 
    barwidth=0.9, title=sprintf("%s putative targets (%s genes)", motif_short, length(pos_top_genes)), enrichLimits=c(0, 6))
  pdf(paste0('./TFregulators', sprintf("/%s_LS.pdf", motif_short)), width=8, height=6)
  print(p)
  print(up_go_plot)
  dev.off()
}

###############内皮细胞Endo######################
Inflamatory <- db %>% filter(gs_name == "HALLMARK_INFLAMMATORY_RESPONSE") %>% pull(human_gene_symbol) 
Apoptosis <- db %>% filter(gs_name == "KEGG_APOPTOSIS") %>% pull(human_gene_symbol)
Chemokine <- db %>% filter(gs_name == "KEGG_CHEMOKINE_SIGNALING_PATHWAY") %>% pull(human_gene_symbol)
TGF <- db %>% filter(gs_name == "HALLMARK_TGF_BETA_SIGNALING") %>% pull(human_gene_symbol)
endo_deg <- subset(deg,celltype=='Endo')
head(endo_deg)
markerGenes <- intersect(c(Inflamatory, Apoptosis,Chemokine,TGF),endo_deg$gene) %>% unique() %>% sort()
markerGenes <- c(markerGenes,"EMCN","PTPRB","FLT1","MEIS2","PECAM1","KDR") %>% unique()
markerGenes
####内皮细胞转录因子
regulators <- c("ETV1_568","ETV4_569","ELF1_565","IKZF1_397","GABPA_576")
# Store results for each TF
res_list <- list()

##########################################################################################
for(motif in regulators){
  motif_short <- strsplit(motif,"_")[[1]][1]
  # First get motif positions
  motifLocs <- motifGR[motifGR$motifName == motif]
  # Calculate Linkage Score for motif
  LS <- calculateLinkageScore(motifLocs, p2gGR)
  # Get just genes correlated to motif
  motifGeneCorDF <- geneCorMat[geneCorMat$motifName == motif,]
  plot_df <- merge(LS, motifGeneCorDF, by="symbol", all.x=TRUE)
  # Calculate motif enrichment per gene
  ME <- calculateMotifEnrichment(motifLocs, p2gGR)
  plot_df <- merge(plot_df, ME, by="symbol", all.x=TRUE)
  plot_df <- plot_df[,c("symbol", "LS", "Correlation", "FDR", "mLog10pval")]
  plot_df$toLabel <- "NO"
  topN <- 5
  plot_df <- plot_df[order(plot_df$LS, decreasing=TRUE),]
  plot_df$rank_LS <- 1:nrow(plot_df)
  plot_df$toLabel[1:topN] <- "YES"
  plot_df <- plot_df[order(plot_df$Correlation, decreasing=TRUE),]
  plot_df$rank_Corr <- 1:nrow(plot_df)
  plot_df$toLabel[1:topN] <- "YES"
  plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=TRUE),]
  plot_df$rank_Pval <- 1:nrow(plot_df)
  plot_df$toLabel[1:10] <- "YES"
  plot_df$meanRank <- apply(plot_df[,c("rank_LS", "rank_Corr", "rank_Pval")], 1, mean)
  plot_df <- plot_df[order(plot_df$meanRank, decreasing=FALSE),]
  plot_df$toLabel[1:topN] <- "YES"
  # Label any marker genes in window of interest
  LS_window <- quantile(plot_df$LS, 0.8)
  corr_window <- 0.25
  pos_top_genes <- plot_df[plot_df$LS > LS_window & plot_df$Correlation > corr_window,]$symbol
  write.csv(pos_top_genes, paste0('./TFregulators', sprintf("/regulatory_targets_%s.csv", motif_short)))
  #neg_top_genes <- plot_df[plot_df$LS > LS_window & -plot_df$Correlation > corr_window,]$symbol
  if(nrow(plot_df[plot_df$symbol %in% pos_top_genes & plot_df$symbol %in% markerGenes,]) > 0){
    plot_df[plot_df$symbol %in% pos_top_genes & plot_df$symbol %in% markerGenes,]$toLabel <- "YES"
  }
  res_list[[motif_short]] <- pos_top_genes # Save regulatory targets
  # Save dataframe of results
  save_df <- plot_df[plot_df$symbol %in% pos_top_genes,c(1:5)]
  save_df <- save_df[order(save_df$Correlation, decreasing=TRUE),]
  saveRDS(save_df, paste0('./TFregulators', sprintf("/regulatory_targets_%s.rds", motif_short)))
  plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=FALSE),]
  # Label motif as well
  plot_df$toLabel[which(plot_df$symbol == motif_short)] <- "YES"
  plot_df$symbol[which(plot_df$toLabel == "NO")] <- ""
  # Threshold pvalue for plotting
  maxPval <- 5
  plot_df$mLog10pval <- ifelse(plot_df$mLog10pval > maxPval, maxPval, plot_df$mLog10pval)
  #Plot results
  p <- (
    ggplot(plot_df, aes(x=Correlation, y=LS, color=mLog10pval)) 
      #+ geom_point(size = 2)
      + ggrastr::geom_point_rast(size=2)
      + ggrepel::geom_text_repel(
          data=plot_df[plot_df$toLabel=="YES",], aes(x=Correlation, y=LS, label=symbol), 
          #data = plot_df, aes(x=Correlation, y=LS, label=symbol), #(don't do this, or the file will still be huge...)
          size=3,
          point.padding=0, # additional pading around each point
          box.padding=0.5,
          min.segment.length=0, # draw all line segments
          max.overlaps=Inf, # draw all labels
          #nudge_x = 2,
          color="black"
      ) 
      + geom_vline(xintercept=0, lty="dashed") 
      + geom_vline(xintercept=corr_window, lty="dashed", color="red")
      + geom_hline(yintercept=LS_window, lty="dashed", color="red")
      + theme_BOR(border=FALSE)
      + theme(panel.grid.major=element_blank(), 
              panel.grid.minor= element_blank(), 
              plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
              aspect.ratio=1.0,
              #legend.position = "none", # Remove legend
              axis.text.x = element_text(angle=90, hjust=1))
      + ylab("Linkage Score") 
      + xlab("Motif Correlation to Gene") 
      + scale_color_gradientn(colors=cmaps_BOR$zissou, limits=c(0, maxPval))
      + scale_y_continuous(expand = expansion(mult=c(0,0.05)))
      + scale_x_continuous(limits = c(-0.85, 0.955)) # Force plot limits
      + ggtitle(sprintf("%s putative targets", motif_short))
      )
  # Positively regulated genes:
  upGO <- calcTopGo(allGenes, interestingGenes=pos_top_genes, nodeSize=5, ontology="BP")
  upGO <- upGO[order(as.numeric(upGO$pvalue), decreasing=FALSE),]
  up_go_plot <- topGObarPlot(upGO, cmap=cmaps_BOR$comet, nterms=3, border_color="black", 
    barwidth=0.9, title=sprintf("%s putative targets (%s genes)", motif_short, length(pos_top_genes)), enrichLimits=c(0, 6))
  pdf(paste0('./TFregulators', sprintf("/%s_LS.pdf", motif_short)), width=8, height=6)
  print(p)
  print(up_go_plot)
  dev.off()
}

###调控网络数据准备
tf_name <- c("HNF4G","HNF4A","HNF1B")
df_lst <- list()
i=1
for(tf in tf_name){
    df <- read.csv(paste0("./TFregulators/regulatory_targets_",tf,".csv"),row.names=1)
    df_lst[[i]] <- data.frame(TF=tf,target=df$x)
    i=i+1
}
data_pt <- do.call(rbind,df_lst)
table(data_pt$TF)
head(data_pt)

tf_name <- c("ETV1","IKZF1","GABPA")
df_lst <- list()
i=1
for(tf in tf_name){
    df <- read.csv(paste0("./TFregulators/regulatory_targets_",tf,".csv"),row.names=1)
    df_lst[[i]] <- data.frame(TF=tf,target=df$x)
    i=i+1
}
data_endo <- do.call(rbind,df_lst)
table(data_endo$TF)
head(data_endo)

####查看target在RNA数据中是否特异性，以及差异
kidney <- readRDS("/data/work/tangyuanchun/graduation/RNA/result/kidney_harmonyfinal.rds")
#####AUCell值
library(AUCell)
#cells_AUC <- AUCell_run(exprMatrix, geneSets)
cells_rankings <- AUCell_buildRankings(kidney@assays$RNA@data)
gene_sets <- list(HNF4G=data_pt[data_pt$TF=='HNF4G',]$target,HNF4A=data_pt[data_pt$TF=='HNF4A',]$target,
                  HNF1B=data_pt[data_pt$TF=='HNF1B',]$target,ETV1=data_endo[data_endo$TF=='ETV1',]$target,
                  IKZF1=data_endo[data_endo$TF=='IKZF1',]$target,GABPA=data_endo[data_endo$TF=='GABPA',]$target)
cells_AUC <- AUCell_calcAUC(gene_sets, cells_rankings, 
                            aucMaxRank=nrow(cells_rankings)*0.1)

#也可以直接使用 AUCell_run 函数得到上述2个步骤相同的结果。
#cells_AUC2 <- AUCell_run(exprMatrix, geneSets)
kidney$HNF4G <- as.numeric(getAUC(cells_AUC)['HNF4G', ])
kidney$HNF4A <- as.numeric(getAUC(cells_AUC)['HNF4A', ])
kidney$HNF1B <- as.numeric(getAUC(cells_AUC)['HNF1B', ])
kidney$ETV1 <- as.numeric(getAUC(cells_AUC)['ETV1', ])
kidney$IKZF1 <- as.numeric(getAUC(cells_AUC)['IKZF1', ])
kidney$GABPA <- as.numeric(getAUC(cells_AUC)['GABPA', ])
plot_df <- kidney@meta.data[,c("HNF4G","HNF4A","HNF1B","ETV1","IKZF1","GABPA","celltype")]
head(plot_df)
#使用aggregate函数取平均
result=aggregate(.~celltype,plot_df,mean)
rownames(result) <- result$celltype
result <- t(result[,-1])
result
#对数据按照基因scale
Result_matrix = t(scale(t(result)))
#Result_matrix = MinMax(Result_matrix, min = -4, max = 4)
levels=c("Immune", "PODO", "ICB", "ICA", "Myofib", "Endo", "DCT", 
         "CNT", "PC", "PEC", "TAL","PT_VCAM1","PT")
Result_matrix <- Result_matrix[,rev(levels)]
library(ComplexHeatmap)
p1 <- Heatmap(Result_matrix,name="Scores",
         col=c("#4264a3","white","#9b1a28"),
         cluster_rows = F,
         cluster_columns = F,
         show_column_names = T,
         show_row_names = T,
         column_title = "TF target gene score",
         row_names_side = "left",column_names_side = "bottom",
         column_names_max_height = unit(200, "cm"))
pdf("./Plots/TF_target_genescore_heatmap.pdf",width = 6,height = 3.5)
p1
dev.off()

##近端小管细胞调控网络
deg1 <- subset(pt_deg,abs(avg_log2FC)>0.5 & p_val_adj<0.05 & (pct.1>0.5 | pct.2>0.5))
dim(deg1)
data1 <- merge(data_pt,deg1,by.x='target',by.y='gene')
head(data1)
write.csv(data1[,c(1,2)],"./TFregulators/PT-tf-target-network.csv",quote=FALSE)

###基因功能富集
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library(ggplot2)
##功能注释
ids=bitr(data1$target,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID

go_all<- compareCluster(ids,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05
)

go_all  = setReadable(go_all, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
go_data<-as.data.frame(go_all)
head(go_data)

metabolic <- unique(unlist(strsplit(go_data[grep("metabolic process",go_data$Description),]$geneID,'/')))
transport <- unique(unlist(strsplit(go_data[grep("transport",go_data$Description),]$geneID,'/')))
intersect(metabolic,transport)

net_lable <- data.frame(key=unique(c(data1$target,data1$TF)))
net_lable$IsTF <- 'NO'
net_lable[net_lable$key %in% unique(data1$TF),]$IsTF <- 'YES'
head(net_lable)
table(net_lable$IsTF)
data2 <- merge(net_lable,deg1,by.x='key',by.y='gene',all.x=TRUE)
head(data2)
dim(data2)
data2[is.na(data2)] <- 0
data2$label <- '1'
data2[data2$avg_log2FC>0 & data2$IsTF=='NO',]$label <- 'up'
data2[data2$avg_log2FC<0 & data2$IsTF=='NO',]$label <- 'down'
table(data2$label)
data2 <- data2[,c(1,2,4,9)]
data2$avg_log2FC[data2$avg_log2FC == 0 | data2$IsTF=='YES'] <- NA
data2$label[data2$label == '1'] <- NA
data2$avg_log2FC <- abs(data2$avg_log2FC)
data2$pathway <- '1'
data2[data2$key %in% metabolic,]$pathway <- "Metabolic"
data2[data2$key %in% transport,]$pathway <- "Transport"
data2$pathway[data2$IsTF=='YES'] <- NA
head(data2)
table(data2$pathway)
write.csv(data2,"./TFregulators/PT-tf-target-network_lable.csv",quote=FALSE)

####内皮细胞调控网络
deg1 <- subset(endo_deg,abs(avg_log2FC)>0.5 & p_val_adj<0.05  & (pct.1>0.3 | pct.2>0.3))
dim(deg1)
data1 <- merge(data_endo,deg1,by.x='target',by.y='gene')
head(data1)
table(data1$TF)
write.csv(data1[,c(1,2)],"./TFregulators/Endo-tf-target-network.csv",quote=FALSE)
###基因功能富集
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library(ggplot2)
##功能注释
ids=bitr(data1$target,'SYMBOL','ENTREZID','org.Hs.eg.db') ## 将SYMBOL转成ENTREZID

go_all<- compareCluster(ids,
                       fun = "enrichGO",
                       OrgDb = "org.Hs.eg.db",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05
)

go_all  = setReadable(go_all, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
go_data<-as.data.frame(go_all)
head(go_data)

proliferation <- unique(unlist(strsplit(go_data[grep("proliferation",go_data$Description),]$geneID,'/')))
chemotaxis <- unique(unlist(strsplit(go_data[grep("chemotaxis",go_data$Description),]$geneID,'/')))
intersect(proliferation,chemotaxis)

net_lable <- data.frame(key=unique(c(data1$target,data1$TF)))
net_lable$IsTF <- 'NO'
net_lable[net_lable$key %in% unique(data1$TF),]$IsTF <- 'YES'
head(net_lable)
table(net_lable$IsTF)
data2 <- merge(net_lable,deg1,by.x='key',by.y='gene',all.x=TRUE)
head(data2)
dim(data2)
data2[is.na(data2)] <- 0
data2$label <- '1'
data2[data2$avg_log2FC>0 & data2$IsTF=='NO',]$label <- 'up'
data2[data2$avg_log2FC<0 & data2$IsTF=='NO',]$label <- 'down'
table(data2$label)
data2 <- data2[,c(1,2,4,9)]
data2$avg_log2FC[data2$avg_log2FC == 0 | data2$IsTF=='YES'] <- NA
data2$label[data2$label == '1'] <- NA
data2$avg_log2FC <- abs(data2$avg_log2FC)
head(data2)
data2$pathway <- '1'
data2[data2$key %in% proliferation,]$pathway <- "Proliferation"
data2[data2$key %in% chemotaxis,]$pathway <- "Chemotaxis"
data2$pathway[data2$IsTF=='YES'] <- NA
head(data2)
table(data2$pathway)
write.csv(data2,"./TFregulators/Endo-tf-target-network_lable.csv",quote=FALSE)