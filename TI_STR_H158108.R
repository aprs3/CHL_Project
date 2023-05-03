setwd("~/Scrivania/CHL_Project")

options(parallelly.fork.enable = TRUE)
library(Matrix)
library(data.table)
library(ggplot2)
library(ggrepel)
library(zeallot)
#devtools::install_github("seriph78/COTAN", ref = "devel")
library(COTAN)

library(rlang)
library(dendextend)
#library(idyr)
library(grid) 
library(ComplexHeatmap)
library(circlize)

source("utils.R")

################################################################################

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
rownames(mat) = feature.names$V2 #one for each feature

mat = as(mat, "dgCMatrix")

obj = COTAN(raw = mat)
obj = initializeMetaDataset(obj, GEO = "", sequencingMethod = "10X", sampleCondition = "")

#stampa il numero di cellule per ogni paziente
clean_cellnames <- c()
for ( col in 1:obj@raw@Dim[[2]]){
  clean_cellname = sub("_.*", "", obj@raw@Dimnames[[2]][[col]])
  #print("_.*", "", obj@raw@Dimnames[[2]][[col]])
  #print(clean_cellname)
  clean_cellnames <-  c(clean_cellnames, clean_cellname)
}

sort(table(clean_cellnames))

patientID = "H158108"

cells_to_remove <- getCells(obj)[!grepl(patientID, colnames(mat))]
#cells_to_remove <- getCells(obj)[!(grepl("I130064|N130064|H158108" , colnames(mat)))]

obj <- dropGenesCells(obj, cells = cells_to_remove)

################################################################################
#plots
ECDPlot(obj, yCut = 300)
cellSizePlot(obj) ##somma del numero di geni contati per ogni cellula (per ogni colonna (cellula), somma i valori righe (gene count)) 
genesSizePlot(obj) ##numeri di geni > 0
mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Percentuale di geni = MT (mitocondriali)
mit[["plot"]]

cells_to_rem <- getCells(obj)[getCellsSize(obj) > 6000]
obj <- dropGenesCells(obj, cells = cells_to_rem)
cellSizePlot(obj) ##somma del numero di geni contati per ogni cellula (per ogni colonna (cellula), somma i valori righe (gene count)) 

columnCounts <- diff(obj@raw@p)
names(columnCounts) <- obj@raw@Dimnames[[2]]
cellGeneNumber <- sort(columnCounts, decreasing = FALSE) 
cells_to_rem <- names(cellGeneNumber)[cellGeneNumber > 2000]
obj <- dropGenesCells(obj, cells = cells_to_rem)
genesSizePlot(obj)

mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Analysis again
to_rem <- mit[["sizes"]][["mit.percentage"]] > 7.5
cells_to_rem <- rownames(mit[["sizes"]])[to_rem]
obj <- dropGenesCells(obj, cells = cells_to_rem)
mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-")
mit[["plot"]]

#rimuove tutte le cellule mitocondriali
genes_to_rem = getGenes(obj)[grep('^MT-', getGenes(obj))] 
cells_to_rem = getCells(obj)[which(getCellsSize(obj) == 0)]
obj = dropGenesCells(obj, genes_to_rem, cells_to_rem)

################################################################################
#Clean iterations
setLoggingLevel(newLevel = 2L) #makes the logging a little more sophisticated
obj <- clean(obj)
cleanPlots <- cleanPlots(obj)

cleanPlots$pcaCells
cleanPlots$genes
cleanPlots$UDE
cleanPlots$nu

##REMOVING CLUSTER B IN DATAFRAME
cells_to_rem <- rownames(cleanPlots$pcaCellsData)[cleanPlots$pcaCellsData[["groups"]] == "B"]
obj <- dropGenesCells(obj, cells = cells_to_rem)

obj <- clean(obj)
cleanPlots <- cleanPlots(obj)

cleanPlots$pcaCells
cleanPlots$genes
cleanPlots$UDE
cleanPlots$nu

UDEPlot = cleanPlots$UDE
nuPlot = cleanPlots$nu

#NB: nu sarebbe ciò che nel paper è descritto come DEI o v_c
nuDf = data.frame("nu" = sort(getNu(obj)), "n" = seq_along(getNu(obj)))
yset = 1 # the threshold to remove low UDE cells
cells_to_rem = rownames(nuDf)[nuDf[["nu"]] < yset]
obj <- dropGenesCells(obj, cells = cells_to_rem)

obj <- clean(obj)
cleanPlots <- cleanPlots(obj)

cleanPlots$pcaCells
cleanPlots$genes
cleanPlots$UDE
cleanPlots$nu

################################################################################
#COTAN ANALYSIS
obj = estimateDispersionBisection(obj, cores = 2)

obj <- calculateCoex(obj)

# saving the structure
t <- dataset_name
outDir <- paste0(getwd(), "/", dataset_name, "/", patientID)
#outDir <- tempdir()
#saveRDS(obj, file = file.path(outDir, paste0(t, ".cotan.RDS")))

################################################################################
quant.p = calculateGDI(obj)
head(quant.p)

#ottenuti dal paper
knownCells = list(
  "Activated fibroblasts CCL19+ ADAMADEC1+" = c('ADAMDEC1', 'CCL19'),
  "Endothelial cells CA4+ CD36+"            = c('CD36', 'CD4'),
  "Endothelial cells CD36+"                 = c('CD36', 'COL15A1', 'PLVAP', 'RBP7', 'TMEM88'),
  "Endothelial cells DARC+"                 = c('C2CD4B', 'CPE', 'DARC', 'GPR126', 'SELE'),
  "Endothelial cells LTC4S+ SEMA3G+"        = c('C10orf10', 'LTC4S', 'SEMA3G'),
  "Fibroblasts ADAMDEC1+"                   = c('ADAMDEC1', 'CCL11', 'CCL13', 'HAPLN1'),
  "Fibroblasts KCNN3+ LY6H+"                = c('C7', 'DPT', 'KCNN3', 'LY6H', 'SCN7A'),
  "Fibroblasts NPY+ SLITRK6+"               = c('EDNRB', 'F3', 'NPY', 'NSG1', 'SLITRK6'),
  "Fibroblasts SFRP2+ SLPI+"                = c('IGFBP6', 'MFAP5', 'SFRP2', 'SLPI'),
  "Fibroblasts SMOC2+ PTGIS+"               = c('ADAMTSL3', 'F3', 'PCSK6', 'PTGIS', 'SMOC2'),
  "Glial cells"                             = c('CDH19', 'GPM6B', 'S100B', 'PLP1', 'NRXN1', 'SCN7A', 'LGI4', 'SPP1'),
  "Lymphatics"                              = c('CCL21', 'MMRN1', 'LYVE1', 'TFPI', 'PPFIBP1'),
  "Myofibroblasts GREM1+ GREM2+"            = c('GREM1', 'GREM2', 'ACTG2', 'DES', 'TAGLN', 'MYH11'),
  "Myofibroblasts HHIP+ NPNT+"              = c('HHIP', 'NPNT', 'SOSTDC1', 'ACTG2', 'ACTA2', 'MYH11', 'TAGLN'),
  "Pericytes HIGD1B+ STEAP4+"               = c('NOTCH3', 'HIGD1B', 'STEAP4', 'COX4I2', 'FABP4'),
  "Pericytes RERGL+ NTRK2+"                 = c('NTRK2', 'RERGL', 'PLN', 'NOTCH3')
)

GDIPlot = GDIPlot(obj, cond = "TI_STR", genes = knownCells)
plot(GDIPlot)
################################################################################

c(gSpace, eigPlot, pcaClustersDF, treePlot) %<-%
  establishGenesClusters(obj, groupMarkers = knownCells,
                         numGenesPerMarker = 25, kCuts = 6)

plot(eigPlot)
plot(treePlot)
UMAPPlot(pcaClustersDF[, 1:10], 
         clusters = pcaClustersDF[["hclust"]],
         elements = knownCells,
         title = "Genes' clusters UMAP Plot")

################################################################################
## Uniform Clustering
'
finds a Uniform clusterizations by means of the GDI. Once a preliminary clusterization
is obtained from the Seurat package methods, each cluster is checked for uniformity
via the function checkClusterUniformity(). Once all clusters are checked, all cells
from the non-uniform clusters are pooled together for another iteration of the
entire process, until all clusters are deemed uniform. In the case only a few cells
are left out (<=50), those are flagged as "not_clustered" and the process is stopped.
'
fineClusters <- cellsUniformClustering(obj, GDIThreshold = 1.4, cores = 2, saveObj = TRUE, outDir = outDir)
obj <- addClusterization(obj, clName = "FineClusters", clusters = fineClusters)

c(coexDF, pValueDF) %<-% DEAOnClusters(obj, clusters = fineClusters)
obj <- addClusterizationCoex(obj, clName = "FineClusters", coexDF = coexDF)

c(mergedClusters, coexDF, pValueDF) %<-%
  mergeUniformCellsClusters(obj, GDIThreshold = 1.4, cores = 2, saveObj = TRUE, outDir = outDir)
obj <- addClusterization(obj, clName = "MergedClusters", clusters = mergedClusters, coexDF = coexDF)

mergedUMAPPlot <- UMAPPlot(coexDF, elements = knownCells, title = "Fine Cluster UMAP Plot")
plot(mergedUMAPPlot)

################################################################################
# Making sense of the clusters
library(tuple)

clustersSummaryPlot(obj, clName = "FineClusters")
clustersSummaryPlot(obj, clName = "MergedClusters")

c(cl, gene, score, pVal, adjPval, DEA, isMarker) %<-% findClustersMarkers(
  obj,
  n = 10,
  clusters = mergedClusters,
  markers = knownCells,
  coexDF = coexDF,
  pValueDF = pValueDF,
  deltaExp = NULL,
  method = "bonferroni"
)

markers <- cl
names(markers) <- gene
knownMarkersList = sort(markers[isMarker == 1])
knownMarkersList

plot1 <- UMAPPlot(t(obj@raw), clusters = mergedClusters, title = "Merged Cluster UMAP Plot")
plot(plot1)

treePlot = clustersTreePlot(obj, clName = "MergedClusters")
plot(treePlot$dend)

c(fullcl, fullgene, fullscore, fullpVal, fulladjPval, fullDEA, fullisMarker) %<-% findClustersMarkers(
  obj,
  n = 5L,
  clusters = mergedClusters,
  coexDF = coexDF,
  pValueDF = pValueDF,
  deltaExp = NULL,
  method = "bonferroni"
)

fullmarkers <- fullcl
names(fullmarkers) <- fullgene
fullmarkers

names(fulladjPval) <- fullgene
fulladjPval
###############################################################################
clustersMarkersHeatmapPlotB(
  obj,
  groupMarkers = knownCells,
  clName = "MergedClusters" 
)[["heatmapPlot"]]



#genesHeatmapPlot(obj, primaryMarkers = knownCells, pValueThreshold = 0.001, symmetric = TRUE)

#gene enrichment
groupMarkers <- c(scan("enrichment_list.txt", what="", sep="\n"))
names(groupMarkers) <- groupMarkers
i = 1
for(element in groupMarkers)
{
  groupMarkers[[i]] <- as.list(groupMarkers[[i]])
  
  i = i + 1
}
groupMarkers <- as.list(groupMarkers)
'
clName <- getClusterizationName(obj, clName = "MergedClusters")
expressionCl <- clustersDeltaExpression(obj, clName = "MergedClusters")
enrichment <- geneSetEnrichment(groupMarkers = groupMarkers, clustersCoex = expressionCl)
'
c(enrichmentHm, enrichmentHmUnclustered, scoreDF) %<-% EnrichmentHeatmap(obj, groupMarkers = groupMarkers, clName = "MergedClusters")
