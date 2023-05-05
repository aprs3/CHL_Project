setwd("~/Scrivania/CHL_Project")

options(parallelly.fork.enable = TRUE)
library(Matrix)
library(data.table)
library(ggplot2)
library(ggrepel)
library(zeallot)
library(COTAN)
library(Seurat)
library(tuple)

source("utils.R")
################################################################################
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

logThis("Loading the gene enrichment data")
groupMarkers <- c(scan("enrichment_list.txt", what="", sep="\n"))
names(groupMarkers) <- groupMarkers
i = 1
for(element in groupMarkers)
{
  groupMarkers[[i]] <- as.list(groupMarkers[[i]])
  
  i = i + 1
}
groupMarkers <- as.list(groupMarkers)

################################################################################
setLoggingLevel(newLevel = 2L) #makes the logging a little more sophisticated

logThis("################################################################################")

args<-commandArgs(TRUE)
args<-c("TI_STR", "N130064", -1, -1, -1) #ATTENZIONE: usare solo per debug

dataset_name = args[1]
dataset_path = paste0(getwd(), "/", dataset_name, "/")

print(paste0("Loading the matrix from path", dataset_path))

mat <- Read10X(data.dir = dataset_path)

logThis("Matrix loaded. Loading the COTANObject")

obj = COTAN(raw = mat)
obj = initializeMetaDataset(obj, GEO = "", sequencingMethod = "10X", sampleCondition = "")

logThis("################################################################################")
#stampa il numero di cellule per ogni paziente
clean_cellnames <- c()
for ( col in 1:obj@raw@Dim[[2]]){
  clean_cellname = sub("_.*", "", obj@raw@Dimnames[[2]][[col]])
  clean_cellnames <-  c(clean_cellnames, clean_cellname)
}

logThis("Cells count per patient: ")
sorted_patiens <- sort(table(clean_cellnames))
print(sorted_patiens[sorted_patiens > 1000])

patientID = args[2]

logThis(paste0("Extracting the cells from patient", patientID))
cells_to_remove <- getCells(obj)[!grepl(patientID, colnames(mat))]
obj <- dropGenesCells(obj, cells = cells_to_remove)

logThis("################################################################################")
logThis("Plotting the raw data")
#Initial plots
ECDPlot(obj, yCut = 300)
cellSizePlot(obj)
genesSizePlot(obj)
mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-")
mit[["plot"]]

logThis("Raw data plotted")

###############################################################################
#cells count outlier removal
cells_count <- args[3]
while(cells_count < 0)
{
  logThis("Please insert a valid value for the cells count threshold")
  
  if(cells_count > 0)
  {
    cells_to_rem <- getCells(obj)[getCellsSize(obj) > cells_count]
    obj <- dropGenesCells(obj, cells = cells_to_rem)
    cellSizePlot(obj) ##somma del numero di geni contati per ogni cellula (per ogni colonna (cellula), somma i valori righe (gene count)) 
    
    answer = readline(prompt="Are you ok with this cut? Do you want to cut further? [y/n]")
    if(answer != 'y')
      cells_count = -1
  }
}

###############################################################################
#Gene count outlier removal
gene_count <- args[4]
while(gene_count < 0)
{
  logThis("Please insert a valid value for the gene count threshold")
  
  if(gene_count > 0)
  {
    columnCounts <- diff(obj@raw@p)
    names(columnCounts) <- obj@raw@Dimnames[[2]]
    cellGeneNumber <- sort(columnCounts, decreasing = FALSE) 
    cells_to_rem <- names(cellGeneNumber)[cellGeneNumber > gene_count]
    obj <- dropGenesCells(obj, cells = cells_to_rem)
    genesSizePlot(obj)
    
    answer = readline(prompt="Are you ok with this cut? Do you want to cut further? [y/n]")
    if(answer != 'y')
      gene_count = -1
  }
}

###############################################################################
#Mitocondrial cells removal
mitocondrial_count <- args[5]
while(mitocondrial_count < 0)
{
  logThis("Please insert a valid value for the mitocondrial threshold")
  
  if(mitocondrial_count > 0)
  {
    mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Analysis again
    to_rem <- mit[["sizes"]][["mit.percentage"]] > mitocondrial_count
    cells_to_rem <- rownames(mit[["sizes"]])[to_rem]
    obj <- dropGenesCells(obj, cells = cells_to_rem)
    mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-")
    mit[["plot"]]
    
    answer = readline(prompt="Are you ok with this cut? Do you want to cut further? [y/n]")
    if(answer != 'y')
      mitocondrial_count = -1
  }
}

logThis("Removing the mitocondrial genes altogether")
genes_to_rem = getGenes(obj)[grep('^MT-', getGenes(obj))] 
cells_to_rem = getCells(obj)[which(getCellsSize(obj) == 0)]
obj = dropGenesCells(obj, genes_to_rem, cells_to_rem)

logThis(paste0("FINAL LOG:"))

logThis("################################################################################")
logThis("Running clean")
obj <- clean(obj)

logThis("Running cleanPlots")
cleanPlots <- cleanPlots(obj)

logThis("Data plotting")
cleanPlots$pcaCells
cleanPlots$genes
cleanPlots$UDE
cleanPlots$nu

###############################################################################
##REMOVING CLUSTER B IN DATAFRAME
remove_clusterB <- args[6]

'if args[5] < 0, it means that the user wants to choose on that specific moment
whether to remove the cluster or not basing on the nu plot
'
if(remove_clusterB < 0)
{
  answer = readline(prompt="Basing on the plots shown, do you want to remove the B cluster? [y/n]")
  if(answer == 'y')
    remove_clusterB = 1
  else
    remove_clusterB = 0
}

if(remove_clusterB == 1)
{
  cells_to_rem <- rownames(cleanPlots$pcaCellsData)[cleanPlots$pcaCellsData[["groups"]] == "B"]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  
  obj <- clean(obj)
  cleanPlots <- cleanPlots(obj)
  
  cleanPlots$pcaCells
  cleanPlots$genes
  cleanPlots$UDE
  cleanPlots$nu
}

###############################################################################
UDEPlot = cleanPlots$UDE
nuPlot = cleanPlots$nu

nu_threshold <- args[7]
while(nu_threshold <= 0)
{
  logThis("Please insert a valid value for the nu threshold")
  
  if(nu_threshold > 0)
  {
    #NB: nu sarebbe ciò che nel paper è descritto come DEI o v_c
    nuDf = data.frame("nu" = sort(getNu(obj)), "n" = seq_along(getNu(obj)))
    cells_to_rem = rownames(nuDf)[nuDf[["nu"]] < nu_threshold]
    obj <- dropGenesCells(obj, cells = cells_to_rem)
    
    obj <- clean(obj)
    cleanPlots <- cleanPlots(obj)
    
    cleanPlots$pcaCells
    cleanPlots$genes
    cleanPlots$UDE
    cleanPlots$nu
    
    answer = readline(prompt="Are you ok with this cut? Do you want to cut further? [y/n]")
    if(answer != 'y')
      nu_threshold = -1
  }
}

logThis("Experiments setup done. The rest is fully automated, so enjoy")

################################################################################
#COTAN ANALYSIS
logThis("Estimating the dispersion bijection")
obj = estimateDispersionBisection(obj, cores = 2)

logThis("Calculating the coex matrices (this might take a while)")
obj <- calculateCoex(obj)

outDir <- paste0(getwd(), "/", dataset_name, "/", patientID)

################################################################################
logThis("GDI calculation")
quant.p = calculateGDI(obj)

GDIPlot = GDIPlot(obj, cond = paste0(dataset_name, ", patient ", patientID), genes = knownCells)
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
         title = paste0("Estabilished genes' clusters UMAP Plot for patient ", dataset_name, ", patient ", patientID))

################################################################################
## Uniform Clustering
logThis("Uniform clustering")
fineClusters <- cellsUniformClustering(obj, GDIThreshold = 1.4, cores = 2, saveObj = TRUE, outDir = outDir)
obj <- addClusterization(obj, clName = "FineClusters", clusters = fineClusters)

logThis("Calculating the fine clusters' coex data")
c(coexDF, pValueDF) %<-% DEAOnClusters(obj, clusters = fineClusters)
obj <- addClusterizationCoex(obj, clName = "FineClusters", coexDF = coexDF)

logThis("Merging the clusters")
c(mergedClusters, coexDF, pValueDF) %<-%
  mergeUniformCellsClusters(obj, GDIThreshold = 1.4, cores = 2, saveObj = TRUE, outDir = outDir)
obj <- addClusterization(obj, clName = "MergedClusters", clusters = mergedClusters, coexDF = coexDF)

logThis("Merging completed. Plotting the data")
mergedUMAPPlot <- UMAPPlot(coexDF, elements = knownCells, title = paste0("Fine Cluster UMAP Plot (", dataset_name, ", patient ", patientID, ")"))
plot(mergedUMAPPlot)

################################################################################
# Making sense of the clusters

clustersSummaryPlot(obj, clName = "FineClusters")
clustersSummaryPlot(obj, clName = "MergedClusters")

logThis("Finding the clusters' cell type using the known cells markers")
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

logThis("The known cell markers found in the data are: ")
markers <- cl
names(markers) <- gene
knownMarkersList = sort(markers[isMarker == 1])
knownMarkersList

logThis("Plotting the UMAP with the clusters informations")
plot1 <- UMAPPlot(t(obj@raw), clusters = mergedClusters, title = paste0("Merged Cluster UMAP Plot (", dataset_name, ", patient ", patientID, ")"))
plot(plot1)

treePlot = clustersTreePlot(obj, clName = title = paste0("Merged Cluster tree Plot (", dataset_name, ", patient ", patientID, ")"))
plot(treePlot$dend)

'
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
'
###############################################################################
logThis("Plotting the cell type heatmap")
clustersMarkersHeatmapPlotB(
  obj,
  groupMarkers = knownCells,
  clName = "MergedClusters" 
)[["heatmapPlot"]]

################################################################################
logThis("Plotting the gene enrichment heatmaps")
c(enrichmentHm, enrichmentHmUnclustered, scoreDF) %<-% EnrichmentHeatmap(obj, groupMarkers = groupMarkers, clName = "MergedClusters")

################################################################################
logThis("Program finished. saving the data.")
logThis("TODO")