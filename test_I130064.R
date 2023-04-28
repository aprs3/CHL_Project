options(parallelly.fork.enable = TRUE)
library(Matrix)
library(data.table)
library(ggplot2)
library(ggrepel)
#devtools::install_github("seriph78/COTAN", ref = "devel")
library(COTAN)

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
#H158108: sano                  messo meglio
#N130064: malato NON infiammato messo meglio
#I130064: malato infiammato     messo meglio
sort(table(clean_cellnames))

patientID = "I130064"

#teniamo solo le cellule del donatore con id N130064 (Chron disease)
toMatch <- c("I130064", "N130064", "H158108")

#cells_to_keep   <- getCells(obj)[grepl("N130064" , colnames(mat))]
cells_to_remove <- getCells(obj)[!(grepl("I130064|N130064|H158108" , colnames(mat)))]

obj <- dropGenesCells(obj, cells = cells_to_remove)


################################################################################
#plots
ECDPlot(obj, yCut = 300)
cellSizePlot(obj) ##somma del numero di geni contati per ogni cellula (per ogni colonna (cellula), somma i valori righe (gene count)) 
genesSizePlot(obj) ##numeri di geni > 0
mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Percentuale di geni = MT (mitocondriali)
mit[["plot"]]

cells_to_rem <- getCells(obj)[getCellsSize(obj) > 4000]
obj <- dropGenesCells(obj, cells = cells_to_rem)
cellSizePlot(obj) ##somma del numero di geni contati per ogni cellula (per ogni colonna (cellula), somma i valori righe (gene count)) 

columnCounts <- diff(obj@raw@p)
names(columnCounts) <- obj@raw@Dimnames[[2]]
cellGeneNumber <- sort(columnCounts, decreasing = FALSE) 
cells_to_rem <- names(cellGeneNumber)[cellGeneNumber > 1500]
obj <- dropGenesCells(obj, cells = cells_to_rem)
genesSizePlot(obj)

mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Analysis again
to_rem <- mit[["sizes"]][["mit.percentage"]] > 10
cells_to_rem <- rownames(mit[["sizes"]])[to_rem]
obj <- dropGenesCells(obj, cells = cells_to_rem)
mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-")
mit[["plot"]]

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
yset = .5 # the threshold to remove low UDE cells
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
t = "TI_STR"
outDir <- tempdir()
saveRDS(obj, file = file.path(outDir, paste0(t, ".cotan.RDS")))

################################################################################
quant.p = calculateGDI(obj)
head(quant.p)

#ottenuti dal paper


knownCells = list(
  "Activated fibroblasts CCL19+ ADAMADEC1+" 	= c('CCL19', 'ADAMDEC1'),
  "Endothelial cells CA4+ CD36+" 	= c('CD4', 'CD36'),
  "Endothelial cells CD36+" 	= c('CD36', 'RBP7', 'TMEM88', 'PLVAP', 'COL15A1'),
  "Endothelial cells DARC+" 	= c('DARC', 'SELE', 'C2CD4B', 'GPR126', 'CPE'),
  "Endothelial cells LTC4S+ SEMA3G+" 	= c('SEMA3G', 'LTC4S', 'C10orf10'),
  "Fibroblasts ADAMDEC1+" 	= c('CCL11', 'ADAMDEC1', 'CCL13', 'HAPLN1'),
  "Fibroblasts KCNN3+ LY6H+" 	= c('KCNN3', 'LY6H', 'DPT', 'C7',  'SCN7A'),
  "Fibroblasts NPY+ SLITRK6+" 	= c('NPY', 'SLITRK6', 'F3', 'EDNRB', 'NSG1'),
  "Fibroblasts SFRP2+ SLPI+" 	= c('SLPI', 'SFRP2', 'IGFBP6', 'MFAP5'),
  "Fibroblasts SMOC2+ PTGIS+" 	= c('SMOC2', 'PTGIS', 'F3', 'PCSK6', 'ADAMTSL3', 'PCSK6'),
  "Glial cells" 	= c('GPM6B', 'S100B', 'PLP1', 'NRXN1', 'CDH19', 'SCN7A', 'LGI4', 'SPP1'),
  "Inflammatory fibroblasts IL11+ CHI3L1+"	= c('CHI3L1', 'IL11', 'MMP3', 'MMP1', 'TNC'),
  "Lymphatics" 	= c('CCL21', 'MMRN1', 'LYVE1', 'TFPI', 'PPFIBP1'),
  "Myofibroblasts GREM1+ GREM2+"	= c('GREM1', 'GREM2', 'ACTG2', 'DES', 'TAGLN', 'MYH11'),
  "Myofibroblasts HHIP+ NPNT+" 	= c('HHIP', 'NPNT', 'SOSTDC1', 'ACTG2', 'ACTA2', 'MYH11', 'TAGLN'),
  "Pericytes HIGD1B+ STEAP4+" 	= c('NOTCH3', 'HIGD1B', 'STEAP4', 'COX4I2', 'FABP4'),
  "Pericytes RERGL+ NTRK2+" 	= c('NTRK2', 'RERGL', 'PLN', 'NOTCH3'),
  "Stromal Cycling cells"	= c('HMGB2', 'UBE2C', 'PTTG1', 'TOP2A', 'MKI67', 'CDC20', 'H2AFZ', 'CCNB1', 'BIRC5', 'NUSAP1')
)

GDIPlot(obj, cond = "TI_STR", genes = knownCells)

################################################################################
'
heatmapPlot(genesLists = knownCells, sets = c(1:3), conditions = c("TI_STR"), dir = outDir)
'

genesHeatmapPlot(obj, primaryMarkers = c('GREM1', 'GREM2', 'ACTG2', 'DES', 'TAGLN', 'MYH11'),
                 pValueThreshold = 0.001, symmetric = TRUE)

estabilishedClusters = establishGenesClusters(obj, groupMarkers = knownCells, numGenesPerMarker = 25, kCuts = 6)
plot(estabilishedClusters$plot.eig)
plot(estabilishedClusters$tree_plot)
plot(estabilishedClusters$g.space)

################################################################################
## Uniform Clustering
'
finds a Uniform clusterizations by means of the GDI. Once a preliminary clusterization
is obtained from the Seurat package methods, each cluster is checked for uniformity
via the function checkClusterUniformity(). Once all clusters are checked, all cells
from the non-uniform clusters are pooled together for another iteration of the
entire process, until all clusters are deemed uniform. In the case only a few cells
are left out (\leq 50≤50), those are flagged as "not_clustered" and the process is stopped.
'
fineClusters <- cellsUniformClustering(obj, GDIThreshold = 1.4, cores = 1, saveObj = TRUE, outDir = outDir)
obj <- addClusterization(obj, clName = "FineClusters", clusters = fineClusters)


DEA <- DEAOnClusters(obj, clusters = fineClusters)
coexDF = DEA$coex
pValueDF = DEA$`p-value`
obj <- addClusterizationCoex(obj, clName = "FineClusters", coexDF = coexDF)


merge <- mergeUniformCellsClusters(obj, GDIThreshold = 1.4, cores = 8, saveObj = TRUE, outDir = outDir)
mergedClusters = merge$clusters
coexDF = merge$coexDF
pValueDF = merge$pValueDF
obj <- addClusterization(obj, clName = "MergedClusters", clusters = mergedClusters, coexDF = coexDF)

mergedUMAPPlot <- UMAPPlot(coexDF, elements = knownCells, title = "Fine Cluster UMAP Plot")
plot(mergedUMAPPlot)
