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
matrix.path <- paste0(matrix_dir, dataset_name, ".scp.matrix.mtx")

mat <- readMM(file = matrix.path)
feature.names = read.delim(features.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)
barcode.names = read.delim(barcode.path,
                           header = FALSE,
                           stringsAsFactors = FALSE)

colnames(mat) = barcode.names$V1 #one for each molecule
rownames(mat) = feature.names$V2 #one for each feature


mat = as(mat, "dgCMatrix")

###CREATING COTAN OBJ###
obj = COTAN(raw = mat)
obj = initializeMetaDataset(obj, GEO = "", sequencingMethod = "10X", sampleCondition = "")

### VIOLIN PLOTS
cellSizePlot(obj) ##somma del numero di geni contati per ogni cellula (per ogni colonna (cellula), somma i valori righe (gene count)) 
genesSizePlot(obj) ##numeri di geni > 0

###MITOCHONDRIAL ANALYSIS
mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Percentuale di geni = MT (mitocondriali)
mit[["plot"]]

##REMOVING CELLS
cells_to_rem <- getCells(obj)[getCellsSize(obj) > 6000]
obj <- dropGenesCells(obj, cells = cells_to_rem)
cellSizePlot(obj) ##somma del numero di geni contati per ogni cellula (per ogni colonna (cellula), somma i valori righe (gene count)) 

##REMOVING GENES
cellGeneNumber <- sort(colSums(as.data.frame(getRawData(obj) > 0)), decreasing = FALSE)
cells_to_rem <- names(cellGeneNumber)[cellGeneNumber > 1250]
obj <- dropGenesCells(obj, cells = cells_to_rem)
gc()
genesSizePlot(obj)

##REMOVING MITOCHONDRIAL
mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Analysis again
to_rem <- mit[["sizes"]][["mit.percentage"]] > 3
cells_to_rem <- rownames(mit[["sizes"]])[to_rem]
obj <- dropGenesCells(obj, cells = cells_to_rem)
mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-")
mit[["plot"]]

mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Analysis again
to_rem <- mit[["sizes"]][["mit.percentage"]] < 1
cells_to_rem <- rownames(mit[["sizes"]])[to_rem]
obj <- dropGenesCells(obj, cells = cells_to_rem)
mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-")
mit[["plot"]]

####DATA CLEANING######
obj <- clean(obj)
cleanPlots <- cleanPlots(obj)

###SHOWING PLOTS
cleanPlots$pcaCellsData
cleanPlots$pcaCells 
cleanPlots$genes
cleanPlots$UDE
cleanPlots$nu

##REMOVING CLUSTER B IN DATAFRAME
cells_to_rem <- rownames(cleanPlots$pcaCellsData)[cleanPlots$pcaCellsData[["groups"]] == "B"]
obj <- dropGenesCells(obj, cells = cells_to_rem)

##CLEANING AGAIN
obj <- clean(obj)
cleanPlots <- cleanPlots(obj)


cleanPlots$pcaCellsData
cleanPlots$pcaCells 
cleanPlots$genes
cleanPlots$UDE
cleanPlots$nu

UDEPlot = cleanPlots$UDE
nuPlot = cleanPlots$nu

nuDf = data.frame("nu" = sort(getNu(obj)), "n" = seq_along(getNu(obj)))
yset = 1.25 # the threshold to remove low UDE cells
cells_to_rem = rownames(nuDf)[nuDf[["nu"]] < yset]
obj <- dropGenesCells(obj, cells = cells_to_rem)

obj <- clean(obj)
cleanPlots <- cleanPlots(obj)

cleanPlots$pcaCellsData
cleanPlots$pcaCells 
cleanPlots$genes
cleanPlots$UDE
cleanPlots$nu



###COTAN ANALYSIS###
obj = estimateDispersionBisection(obj, cores = 10)

# COEX evaluation and storing
obj <- calculateCoex(obj)

# saving the structure
t = "TI_STR"
outDir <- "/home/aprs3/Scrivania/CHL_Project/plots"
saveRDS(obj, file = file.path(outDir, paste0(t, ".cotan.RDS")))

##blocked for genes list


### Uniform Clustering

fineClusters <- cellsUniformClustering(obj, GDIThreshold = 1.4, cores = 10,
                                       saveObj = TRUE, outDir = outDir)
obj <- addClusterization(obj, clName = "FineClusters", clusters = fineClusters)

DEAOnCluster = DEAOnClusters(obj, clusters = fineClusters)
coexDF = DEAOnCluster$coex
pValueDF = DEAOnCluster$`p-value`

obj <- addClusterizationCoex(obj, clName = "FineClusters", coexDF = coexDF)

fineUMAPPlot <- UMAPPlot(coexDF, title = "Fine Cluster UMAP Plot")
plot(fineUMAPPlot)


merged = mergeUniformCellsClusters(obj, GDIThreshold = 1.4, cores = 10, saveObj = TRUE, outDir = outDir)

mergedClusters = merged$clusters
coexDF = merged$coexDF
pValueDF = merged$pValueDF

obj <- addClusterization(obj, clName = "MergedClusters", clusters = mergedClusters, coexDF = coexDF)

mergedUMAPPlot <- UMAPPlot(coexDF, title = "Fine Cluster UMAP Plot")
plot(mergedUMAPPlot)

