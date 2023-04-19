library(Matrix)
library(COTAN)
library(data.table)
library(ggplot2)
library(ggrepel)
setwd("~/Scrivania/CHL_Project")

dataset_name = "TI_STR"
dataset_path = paste0(getwd(), "/", dataset_name, "/")

matrix_dir <- dataset_path
barcode.path <- paste0(matrix_dir, dataset_name, ".scp.barcodes.tsv")
features.path <- paste0(matrix_dir, dataset_name, ".scp.features.tsv")
matrix.path <- paste0(matrix_dir, dataset_name, ".scp.raw.mtx")

mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)

colnames(mat) = barcode.names$V1 #one for each molecule
rownames(mat) = feature.names$V1 #one for each feature

cellls_to_keep=2000

mat2 = t(head(t(mat), cellls_to_keep))
dataframe = as.data.frame(as.matrix(mat2)) #Column molecules, rows genes

colnames(dataframe) = barcode.names$V1[1:cellls_to_keep]
rownames(dataframe) = feature.names$V2


obj = new("scCOTAN", raw = dataframe, n_cells=cellls_to_keep)
obj = initializeMetaDataset(obj, GEO="", sequencingMethod = "10X", sampleCondition = "")



#Removing mithocondrial and cells
#genes_to_rem = rownames(obj@raw[grep("^mt", rownames(obj@raw)), ])
#obj@raw = obj@raw[!rownames(obj@raw) %in% genes_to_rem, ]
#cells_to_rem = colnames(obj@raw[which(colSums(obj@raw) == 0)])
#obj@raw = obj@raw[, !colnames(obj@raw) %in% cells_to_rem]

####DATA CLEANING######
ttm = clean(obj)

plots = cleanPlots(ttm)

pcaCells = plots$pcaCells
genes = plots$genes
UDE = plots$UDE
nu = plots$nu

cellsUniformClustering =cellsUniformClustering(ttm, cores = 14, maxIterations = 1)
seuratClustering = seuratClustering(ttm)


