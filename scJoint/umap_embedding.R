#########################################################################
# File Name: umap_embedding.R
# scjoint结果处理
#########################################################################
# load library
library(ggplot2)
#library(ggthemes)
library(scattermore)
library(uwot)
library(grDevices)
library(ArchR)
library(dplyr)
library(tibble)

args = commandArgs(T)
outdir <- args[1]
setwd(outdir)

# Reading label reference
label_class <- read.delim("./kidney_label_to_idx.txt",header = FALSE)       ###read.delim不要求所有列都对等，会按最大列，或指定的列数填充
label_class_num <- unlist(lapply(strsplit(label_class$V1, " "),
                                 function(x) x[length(x)]))
label_class_name <- unlist(lapply(strsplit(label_class$V1, " "),
                                  function(x) paste(x[-length(x)], collapse = " ")))
label_class <- data.frame(name = label_class_name,num = label_class_num)



# The folder with output
results_dir <- "../output"
embedding_files <- list.files(results_dir, "embeddings.txt")

embedding <- list()
for (i in 1:length(embedding_files)) {
    embedding[[i]] <- read.delim(file.path(results_dir, embedding_files[i]),
                                 header = FALSE, sep = " ")
}

names(embedding) <- gsub("_embeddings.txt", "", embedding_files)

cat("Dimension of embedding: ")
print(lapply(embedding, dim))


# Reading KNN prediction
knn_prediction_files <- list.files(results_dir, pattern = "knn_predictions.txt")

knn_prediction <- list()
for (i in 1:length(knn_prediction_files)) {
    knn_prediction[[i]] <- read.delim(file.path(results_dir, knn_prediction_files[i]),
                                      header = FALSE, sep = " ")
    knn_prediction[[i]] <- label_class$name[knn_prediction[[i]]$V1 + 1]
}

names(knn_prediction) <- gsub("_knn_predictions.txt", "", knn_prediction_files)


rna_dataset <- setdiff(names(embedding), names(knn_prediction))
print(rna_dataset)
rna_prediction <- list()
for (i in 1:length(rna_dataset)) {
    rna_prediction[[i]] <- read.delim(file.path(results_dir, paste0(rna_dataset[i], "_predictions.txt")),
                                      header = FALSE, sep = " ")
    rna_prediction[[i]] <- label_class$name[apply(rna_prediction[[i]], 1, which.max)]
}

names(rna_prediction) <- rna_dataset


prediction_list <- append(rna_prediction, knn_prediction)
prediction_list <- prediction_list[names(embedding)]

batch <- rep(names(prediction_list), unlist(lapply(prediction_list, length)))
combine_embedding <- do.call(rbind, embedding)
prediction <- do.call(c, prediction_list)


idx <- sort(sample(length(batch), round(length(batch))))
combine_embedding <- combine_embedding[idx, ]
prediction <- prediction[idx]
batch <- batch[idx]

cat("Dimension to be visualised: ")
print(dim(combine_embedding))

set.seed(1234)
umap_res <- uwot::umap(combine_embedding,metric = "cosine",nn_method="annoy",
                       n_neighbors = 15,min_dist = 0.2)

df <- data.frame(UMAP1 = umap_res[, 1], UMAP2 = umap_res[, 2],
				 celltype = prediction,
				 batch = batch)

rna <- read.csv("kidney_rna_cell.csv")
atac <- read.csv("kidney_atac_cell.csv")
cell <- c(atac$x,rna$x)
rownames(combine_embedding) <- cell
write.csv(combine_embedding,"combine_embedding.csv")
rownames(df) <- cell

df[grep("atac$",df$batch),]$batch <- "ATAC"
df[grep("rna$",df$batch),]$batch <- "RNA"
df$state='1'
df[grep("^M[1,2]",rownames(df)),]$state <- "CTRL"
df[grep("^M[7,9]",rownames(df)),]$state <- "DM"

write.csv(df, "kidney_umap_embedding.csv")

