setwd("~/Scrivania/CHL_Project")

options(parallelly.fork.enable = TRUE)
library(Matrix)
library(data.table)
library(ggplot2)
library(ggrepel)
library(zeallot)
library(COTAN)
library(Seurat)

source("utils.R")
################################################################################
logThis("################################################################################")
setLoggingLevel(newLevel = 2L) #makes the logging a little more sophisticated

args<-commandArgs(TRUE)
args<-c("TI_IMM", "I139892", -1, -1, -1, -1, -1, -1, -1) #ATTENZIONE: usare solo per debug

dataset_name = args[1]
dataset_path = paste0(getwd(), "/", dataset_name, "/dataset/")

print(paste0("Loading the matrix from path", dataset_path))

mat <- Read10X(data.dir = dataset_path)

logThis("Matrix loaded. Loading the COTANObject")

obj = COTAN(raw = mat)
obj = initializeMetaDataset(obj, GEO = "", sequencingMethod = "10X", sampleCondition = "")

logThis("################################################################################")

#ottenuti dal paper
knownCells = read_csv_data(paste0(getwd(), "/", dataset_name, "/known_cells_genes.csv"))
knownCellsClean <- knownCells

logThis("Loading the gene enrichment data")

###APPENDING GENE ENRICHMENT LIST TO KNOWN CELLS  
groupMarkers <- c(scan("enrichment_list.txt", what="", sep="\n"))
temp <- list(groupMarkers)
names(temp) <- "Enrichment genes"
knownCellsTmp <- append(knownCells,temp)

names(groupMarkers) <- groupMarkers
i = 1
for(element in groupMarkers)
{
  groupMarkers[[i]] <- as.list(groupMarkers[[i]])
  
  i = i + 1
}

groupMarkers <- as.list(groupMarkers)
groupMarkers <- c(knownCells, groupMarkers) ##Genes family (from 1 to n_cells) + list with enrichment genes (332 genes)

knownCells <- knownCellsTmp #Genes family (n_cells) + enrichment genes in one list (1 that contains 332 genes) 
################################################################################
#stampa il numero di cellule per ogni paziente
clean_cellnames <- c()
for ( col in 1:obj@raw@Dim[[2]]){
  clean_cellname = sub("_.*", "", obj@raw@Dimnames[[2]][[col]])
  clean_cellnames <-  c(clean_cellnames, clean_cellname)
}

###PRINTING PATIENTS WITH MORE THAN 1000 CELLS
logThis("Cells count per patient: ")
sorted_patiens <- sort(table(clean_cellnames))
print(sorted_patiens[sorted_patiens > 1000])

patientID = args[2]
outDir <- paste0(getwd(), "/", dataset_name, "/", patientID, "/")
outDirPlot <- paste0(outDir, "/plot/")
outDirClustering <- paste0(outDir, "/clustering/")

#Creating patient folder
dir.create(file.path(paste0(getwd(), "/", dataset_name, "/", patientID)))

#Creating plot folder
dir.create(file.path(paste0(getwd(), "/", dataset_name, "/", patientID, "/plot")))

#Creating clustering folder
dir.create(file.path(paste0(getwd(), "/", dataset_name, "/", patientID, "/clustering")))

logThis(paste0("Extracting the cells from patient", patientID))
cells_to_remove <- getCells(obj)[!grepl(patientID, colnames(mat))] #Removing all patients not considered
obj <- dropGenesCells(obj, cells = cells_to_remove)

rm(mat)
gc()

logThis("################################################################################")
logThis("Plotting the raw data")
#Initial plots
plotPDF(outDirPlot,dataset_name,patientID,"_00_ECDPlot",ECDPlot(obj, yCut = 300),width = NULL, height = NULL)
plotPDF(outDirPlot,dataset_name,patientID,"_01_CellSizePlot",cellSizePlot(obj))
plotPDF(outDirPlot,dataset_name,patientID,"_02_GenesSizePlot",genesSizePlot(obj))
mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-")
plotPDF(outDirPlot,dataset_name,patientID,"_03_MitocondrialPlot",mit[["plot"]])

logThis("Raw data plotted")
###############################################################################
#cells count outlier removal
i = 1 #counter for plots
cells_count <- as.integer(args[3])

if(cells_count >= 0)
{
  cells_to_rem <- getCells(obj)[getCellsSize(obj) > cells_count]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_04_CellSizePlot_cut_", i), cellSizePlot(obj)) ##somma del numero di geni contati per ogni cellula (per ogni colonna (cellula), somma i valori righe (gene count)) 
} else {
  while(cells_count < 0)
  {
    cells_count = as.integer(readline(prompt="Please insert a valid value for the cells count threshold:"))
    
    if(cells_count > 0)
    {
      cells_to_rem <- getCells(obj)[getCellsSize(obj) > cells_count]
      obj <- dropGenesCells(obj, cells = cells_to_rem)
      plotPDF(outDirPlot,dataset_name,patientID,paste0("_04_CellSizePlot_cut_", i),cellSizePlot(obj)) ##somma del numero di geni contati per ogni cellula (per ogni colonna (cellula), somma i valori righe (gene count)) 
      
      answer = readline(prompt="Do you want to cut further? [y/n] ")
      if(answer == 'y')
        cells_count = -1
      
      i = i+1
    }
  }
}


###############################################################################
#Gene count outlier removal
gene_count <- as.integer(args[4])
i = 1
if(gene_count >= 0)
{
  columnCounts <- diff(obj@raw@p)
  names(columnCounts) <- obj@raw@Dimnames[[2]]
  cellGeneNumber <- sort(columnCounts, decreasing = FALSE) 
  cells_to_rem <- names(cellGeneNumber)[cellGeneNumber > gene_count]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_05_GeneCountPlot_cut_", i),genesSizePlot(obj))
} else{
  while(gene_count < 0)
  {
    gene_count = as.integer(readline(prompt="Please insert a valid value for the genes count threshold: "))
    
    if(gene_count > 0)
    {
      columnCounts <- diff(obj@raw@p)
      names(columnCounts) <- obj@raw@Dimnames[[2]]
      cellGeneNumber <- sort(columnCounts, decreasing = FALSE) 
      cells_to_rem <- names(cellGeneNumber)[cellGeneNumber > gene_count]
      obj <- dropGenesCells(obj, cells = cells_to_rem)
      plotPDF(outDirPlot,dataset_name,patientID,paste0("_05_GeneCountPlot_cut_", i),genesSizePlot(obj))
      
      answer = readline(prompt="Do you want to cut further? [y/n] ")
      if(answer == 'y')
        gene_count = -1
      
      i = i+1
    }
  }
}


###############################################################################
#Mitocondrial cells removal
mitocondrial_count <- as.numeric(args[5])
i = 1

if(mitocondrial_count >= 0)
{
  mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Analysis again
  to_rem <- mit[["sizes"]][["mit.percentage"]] > mitocondrial_count
  cells_to_rem <- rownames(mit[["sizes"]])[to_rem]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-")
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_06_MitocondrialCount_cut_", i),mit[["plot"]])
} else {
  while(mitocondrial_count < 0)
  {
    mitocondrial_count = as.numeric(readline(prompt="Please insert a valid value for the mitocondrial percentage threshold: "))
    
    if(mitocondrial_count > 0)
    {
      mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Analysis again
      to_rem <- mit[["sizes"]][["mit.percentage"]] > mitocondrial_count
      cells_to_rem <- rownames(mit[["sizes"]])[to_rem]
      obj <- dropGenesCells(obj, cells = cells_to_rem)
      mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-")
      plotPDF(outDirPlot,dataset_name,patientID,paste0("_06_MitocondrialCount_cut_", i),mit[["plot"]])
      
      answer = readline(prompt="Do you want to cut further? [y/n] ")
      if(answer == 'y')
        mitocondrial_count = -1
      
      i = i+1
    }
  }
}

##REMOVING ALL MT MITOCONDRIAL CELLS
logThis("Removing the mitocondrial genes altogether")
genes_to_rem = getGenes(obj)[grep('^MT-', getGenes(obj))] 
cells_to_rem = getCells(obj)[which(getCellsSize(obj) == 0)]
obj = dropGenesCells(obj, genes_to_rem, cells_to_rem)

logThis("################################################################################")
logThis("Running clean")
obj <- clean(obj)

logThis("Running cleanPlots")
cleanPlots <- cleanPlots(obj)

logThis("Data plotting")
plotPDF(outDirPlot,dataset_name,patientID,("_07_PCACells"),cleanPlots$pcaCells)
###############################################################################
##REMOVING CLUSTER B IN DATAFRAME
remove_clusterB <- as.integer(args[6])

#if args[5] < 0, it means that the user wants to choose on that specific moment
#whether to remove the cluster or not basing on the nu plot
if(remove_clusterB < 0)
{
  answer = readline(prompt="Basing on the plots shown, do you want to remove the B cluster? [y/n] ")
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
  
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_08_PCACellsBRemoval"),cleanPlots$pcaCells)
}

plotPDF(outDirPlot,dataset_name,patientID,"_09_CleanPlotGenes",cleanPlots$genes)
plotPDF(outDirPlot,dataset_name,patientID,"_10_CleanPlotUDE",cleanPlots$UDE)
plotPDF(outDirPlot,dataset_name,patientID,"_11_CleanPlotNu",cleanPlots$nu)
###############################################################################
####NU CLEANING
UDEPlot = cleanPlots$UDE
nuPlot = cleanPlots$nu

nu_threshold <- as.numeric(args[7])
i = 1

if(nu_threshold > 0)
{
  #NB: nu sarebbe ciò che nel paper è descritto come DEI o v_c
  nuDf = data.frame("nu" = sort(getNu(obj)), "n" = seq_along(getNu(obj)))
  cells_to_rem = rownames(nuDf)[nuDf[["nu"]] < nu_threshold]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  
  obj <- clean(obj)
  cleanPlots <- cleanPlots(obj)
  
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_12_PCACells_cut_", i),cleanPlots$pcaCells)
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_13_CleanPlotGenes_cut_", i),cleanPlots$genes)
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_14_CleanPlotUDE_cut_", i),cleanPlots$UDE)
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_15_CleanPlotNu_cut_", i),cleanPlots$nu)
} else {
  while(nu_threshold <= 0)
  {
    nu_threshold = as.numeric(readline(prompt="Please insert a valid value for the nu threshold (insert 0.0000000000000001 for skipping): "))
    
    if(nu_threshold > 0)
    {
      #NB: nu sarebbe ciò che nel paper è descritto come DEI o v_c
      nuDf = data.frame("nu" = sort(getNu(obj)), "n" = seq_along(getNu(obj)))
      cells_to_rem = rownames(nuDf)[nuDf[["nu"]] < nu_threshold]
      obj <- dropGenesCells(obj, cells = cells_to_rem)
      
      obj <- clean(obj)
      cleanPlots <- cleanPlots(obj)
      
      plotPDF(outDirPlot,dataset_name,patientID,paste0("_12_PCACells_cut_", i),cleanPlots$pcaCells)
      plotPDF(outDirPlot,dataset_name,patientID,paste0("_13_CleanPlotGenes_cut_", i),cleanPlots$genes)
      plotPDF(outDirPlot,dataset_name,patientID,paste0("_14_CleanPlotUDE_cut_", i),cleanPlots$UDE)
      plotPDF(outDirPlot,dataset_name,patientID,paste0("_15_CleanPlotNu_cut_", i),cleanPlots$nu)
      
      answer = readline(prompt="Do you want to cut further? [y/n] ")
      if(answer == 'y')
        nu_threshold = -1
      
      i = i + 1
    }
  }
}

################################################################################
#COTAN ANALYSIS
logThis("Estimating the dispersion bijection")
obj = estimateDispersionBisection(obj, cores = 4)

logThis("Calculating the coex matrices (this might take a while)")
obj <- calculateCoex(obj)
################################################################################
logThis("GDI calculation")
quant.p = calculateGDI(obj)

GDIPlot = GDIPlot(obj, cond = paste0("GDI plot for dataset ", dataset_name, ", patient ", patientID), genes = knownCellsClean)
plotPDF(outDirPlot,dataset_name,patientID,"_16_GDIPlot",GDIPlot, width = 14, height = 14)
################################################################################
## Uniform Clustering
logThis("Uniform clustering")
fineClusters <- cellsUniformClustering(obj, GDIThreshold = 1.4, cores = 2, saveObj = TRUE, outDir = outDirClustering)
obj <- addClusterization(obj, clName = "FineClusters", clusters = fineClusters)

logThis("Calculating the fine clusters' coex data")
c(coexDF, pValueDF) %<-% DEAOnClusters(obj, clusters = fineClusters)
obj <- addClusterizationCoex(obj, clName = "FineClusters", coexDF = coexDF)

logThis("Merging the clusters")
c(mergedClusters, coexDF, pValueDF) %<-%
  mergeUniformCellsClusters(obj, GDIThreshold = 1.4, cores = 6, saveObj = TRUE, outDir = outDirClustering)
obj <- addClusterization(obj, clName = "MergedClusters", clusters = mergedClusters, coexDF = coexDF)
################################################################################
# Making sense of the clusters
plotPDF(outDirPlot,dataset_name,patientID,"_22_FineClustersSummary",clustersSummaryPlot(obj, clName = "FineClusters")[["plot"]])
plotPDF(outDirPlot,dataset_name,patientID,"_23_MergedClustersSummary",clustersSummaryPlot(obj, clName = "MergedClusters")[["plot"]])

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

###############################################################################
logThis("Plotting the cell type heatmap")
kCuts <- as.numeric(args[8])

if(kCuts >= 1)
{
  clustersMarkersHeatmapPlot <- clustersMarkersHeatmapPlotB(
    obj,
    kCuts = kCuts,
    groupMarkers = knownCells,
    clName = "MergedClusters")[["heatmapPlot"]]
  plotPDF(outDirPlot,dataset_name,patientID,"_26_clustersMarkersHeatmapPlot",clustersMarkersHeatmapPlot, width=20, height=20)
} else {
  clustersMarkersHeatmapPlot <- clustersMarkersHeatmapPlotB(
    obj,
    kCuts = 1,
    groupMarkers = knownCells,
    clName = "MergedClusters")[["heatmapPlot"]]
  plotPDF(outDirPlot,dataset_name,patientID,"_26_clustersMarkersHeatmapPlot",clustersMarkersHeatmapPlot, width=20, height=20)
  
  while(kCuts < 1)
  {
    kCuts = as.numeric(readline(prompt="Please insert a valid value for the number of cuts: "))
    
    if(kCuts >= 1)
    {
      clustersMarkersHeatmapPlot <- clustersMarkersHeatmapPlotB(
        obj,
        kCuts = kCuts,
        groupMarkers = knownCells,
        clName = "MergedClusters")[["heatmapPlot"]]
      plotPDF(outDirPlot,dataset_name,patientID,"_26_clustersMarkersHeatmapPlot",clustersMarkersHeatmapPlot, width=20, height=20)
      
      answer = readline(prompt="Do you want to cut differently? [y/n]")
      if(answer == 'y')
        kCuts = -1
    }
  }
}

################################################################################
logThis("Plotting the gene enrichment heatmaps")
enrichment_cut <- as.integer(args[9])

if(enrichment_cut >= 1)
{
  c(enrichmentHm, enrichmentHmUnclustered, scoreDF, enrichmentClusters) %<-% EnrichmentHeatmap(obj, row_km = kCuts,  column_km = enrichment_cut, groupMarkers = groupMarkers, clName = "MergedClusters")
  names(enrichmentClusters) = as.character(seq(1, length(enrichmentClusters), by=1))
  save_list_to_csv(enrichmentClusters, outDir, "EnrichmentGenesClusters.csv", header = "cluster_ID,genes")
  
  plotPDF(outDirPlot,dataset_name,patientID,"_27_enrichmentHm",enrichmentHm, width = 35, height = 7.5)
} else {
  c(enrichmentHm, enrichmentHmUnclustered, scoreDF, enrichmentClusters) %<-% EnrichmentHeatmap(obj, row_km = kCuts,  column_km = 1, groupMarkers = groupMarkers, clName = "MergedClusters")
  plotPDF(outDirPlot,dataset_name,patientID,"_27_enrichmentHm",enrichmentHm, width = 35, height = 7.5)
  
  while(enrichment_cut < 1)
  {
    enrichment_cut = as.numeric(readline(prompt="Please insert a valid value for the number of cuts: "))
    
    if(enrichment_cut >= 1)
    {
      c(enrichmentHm, enrichmentHmUnclustered, scoreDF, enrichmentClusters) %<-% EnrichmentHeatmap(obj, row_km = kCuts,  column_km = enrichment_cut, groupMarkers = groupMarkers, clName = "MergedClusters")
      names(enrichmentClusters) = as.character(seq(1, length(enrichmentClusters), by=1))
      save_list_to_csv(enrichmentClusters, outDir, "EnrichmentGenesClusters.csv", header = "cluster_ID,genes")
      
      plotPDF(outDirPlot,dataset_name,patientID,"_27_enrichmentHm",enrichmentHm, width = 35, height = 7.5)
      
      answer = readline(prompt="Do you want to cut differently? [y/n]")
      if(answer == 'y')
        enrichment_cut = -1
    }
  }
}

plotPDF(outDirPlot,dataset_name,patientID,"_28_enrichmentHmUnclustered",enrichmentHmUnclustered, width = 35, height = 5)
################################################################################
logThis("Saving various clusters alternatives")
for(i in 5:20)
{
  c(a, b, c, curr_cluster) %<-% EnrichmentHeatmap(obj, row_km = kCuts,  column_km = i, groupMarkers = groupMarkers, clName = "MergedClusters")
  names(curr_cluster) = as.character(seq(1, length(curr_cluster), by=1))
  save_list_to_csv(curr_cluster, outDir, paste0(dataset_name,"_",patientID,"_","EnrichmentGenesClusters_", i, ".csv"), header = "cluster_ID,genes")
}


################################################################################
logThis(paste0("FINAL LOG:"))
logThis(paste0("cells count: ", cells_count))
logThis(paste0("gene count: ", gene_count))
logThis(paste0("mitocondrial percentage: ", mitocondrial_count))
logThis(paste0("Cluster B removed: ", remove_clusterB))
logThis(paste0("nu_threshold: ", nu_threshold))
logThis(paste0("kCuts: ", kCuts))
logThis(paste0("enrichment_cut ", enrichment_cut))

logThis("Program finished. saving the data.")
save.image(paste0(getwd(), "/", dataset_name, "/", patientID, "/", patientID, "data.RData"))