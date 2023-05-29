#' Performs the data cleaning, clustering, and gene set enrichment analysis
#' using COTAN. Outputs multiple cluster cuts from a given range as .csv files.
#' It can be run by either setting the arguments by command line or by setting
#' the "args" array manually such that the first and second element of
#' the list are respectively the dataset folder name and the patient's ID.
#' The script will create all the folders automatically. If one of the args is
#' set to -1, the script will procedurally ask the user to insert variable
#' values in an interactive way.
#' Also, this script will calculate the clusterization for the targeted enrichment
#' genes for all the number of clusters ranging between two variables, "start" and "end".

#'set it to whatever you like. The working directory path should contain at least
#'the main.R script
setwd("~/Scrivania/CHL_Project")

#'imports the needed libraries
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
#'makes the logging a little more detailed
setLoggingLevel(newLevel = 2L) 

#'Extracts the arguments. In order:
#'1)Name of the folder containing the dataset (it should be inside the working directory)
#'2)ID of the patient
#'3)Cell count threshold. Any cell with higher counts will be removed (-1: insert manually at runtime)
#'4)Gene count threshold. Any cell with higher counts will be removed (-1: insert manually at runtime)
#'5)Mitochondrial percentage threshold. Any cell with higher counts will be removed (-1: insert manually at runtime)
#'6)Whether to remove the "B" cluster (1: yes, anything else: insert manually at runtime)
#'7)Nu threshold. Any cell with lower counts will be removed (-1: insert manually at runtime)
#'8)kCuts argument for the clustersTreePlot method (-1: insert manually at runtime)
#'9)column_km argument for the gene enrichment clusterization (-1: insert manually at runtime)
args<-commandArgs(TRUE)

#alternatively, we can set the args array manually
#args<-c("CO_IMM", "I130084", -1, -1, -1, -1, -1, -1, -1)

#'extracts the dataset name and computes the path toward the dataset folder
#'where the .mtx file is stored (for example, TI_IMM/dataset/)
dataset_name = args[1]
dataset_path = paste0(getwd(), "/", dataset_name, "/dataset/")

print(paste0("Loading the matrix from path", dataset_path))

#'loads the matrix using Seurat (10x file format only, compressed as a .gz folder)
mat <- Read10X(data.dir = dataset_path)

logThis("Matrix loaded. Loading the COTANObject")

#'Creates the COTAN object and loads the data into it
obj = COTAN(raw = mat)
obj = initializeMetaDataset(obj, GEO = "", sequencingMethod = "10X", sampleCondition = "")

logThis("################################################################################")

#'loads the list of known cell types alongside their associated gene markers
knownCells = read_csv_data(paste0(getwd(), "/", dataset_name, "/known_cells_genes.csv"))

#'used only in the GDI plot
knownCellsClean <- knownCells

logThis("Loading the gene enrichment data")

#'loads the list of target genes for the enrichment analysis
groupMarkers <- c(scan("enrichment_list.txt", what="", sep="\n"))

#'creates knownCellsTmp, which is used to calculate the overall enrichment of
#'the target genes in the various clusters (displayed in the last column of the 
#'various "clustersMarkersHeatmapPlot.pdf" files)
temp <- list(groupMarkers)
names(temp) <- "Enrichment genes"
knownCellsTmp <- append(knownCells,temp)

#'casts the various lists as proper lists (otherwise it doesn't work well)
names(groupMarkers) <- groupMarkers
i = 1
for(element in groupMarkers)
{
  groupMarkers[[i]] <- as.list(groupMarkers[[i]])
  
  i = i + 1
}

#'Genes family (from 1 to n_cells) + list with enrichment genes (332 genes)
#'used for the gene set enrichment analysis
groupMarkers <- as.list(groupMarkers)
groupMarkers <- c(knownCells, groupMarkers) 

#'Genes family (n_cells) + enrichment genes in one list (1 that contains 332 genes) 
#'used for the clusterization
knownCells <- knownCellsTmp 
################################################################################
#'Calculate how many cells the dataset have for each patient.
#'Useful for the initial analysis of the dataset. 
clean_cellnames <- c()
for ( col in 1:obj@raw@Dim[[2]]){
  clean_cellname = sub("_.*", "", obj@raw@Dimnames[[2]][[col]])
  clean_cellnames <-  c(clean_cellnames, clean_cellname)
}

#'Prints only the patients with more than 1000 genes
logThis("Cells count per patient: ")
sorted_patiens <- sort(table(clean_cellnames))
print(sorted_patiens[sorted_patiens > 1000])

#'Calculates and generates the paths where we have to store the various analysis results
patientID = args[2]
outDir <- paste0(getwd(), "/", dataset_name, "/", patientID, "/") #The patient's folder inside the dataset's
outDirPlot <- paste0(outDir, "/plot/")                            #The patient plot's folder inside the dataset's
outDirClustering <- paste0(outDir, "/clustering/")                #The patient cluster data's folder inside the dataset's
outDirEnrichmentCsv <- paste0(getwd(), "/", dataset_name, "/", patientID, "/enrichment_csv")
#'Creates the patient's folder
dir.create(file.path(paste0(getwd(), "/", dataset_name, "/", patientID)))

#'Creates the patient plot's folder
dir.create(file.path(paste0(getwd(), "/", dataset_name, "/", patientID, "/plot")))

#'Creates the patient data's folder where we store all the clustering pdfs
dir.create(file.path(paste0(getwd(), "/", dataset_name, "/", patientID, "/clustering")))

#'Creates the patient's folder where we store the .csv(s) of the clusterization of the target genes for the enrichment
dir.create(file.path(outDirEnrichmentCsv))

#'Removes the cells coming from patients other than the one with ID "patientID"
logThis(paste0("Extracting the cells from patient", patientID))
cells_to_remove <- getCells(obj)[!grepl(patientID, colnames(mat))]
obj <- dropGenesCells(obj, cells = cells_to_remove)

#'Frees some memory
rm(mat)
gc()

logThis("################################################################################")
logThis("Plotting the raw data")
#'Plots: The ECD plot, the cell size plot, the genes size plot, and the mitochondrial percentages plot
plotPDF(outDirPlot,dataset_name,patientID,"_00_ECDPlot",ECDPlot(obj, yCut = 300),width = NULL, height = NULL)
plotPDF(outDirPlot,dataset_name,patientID,"_01_CellSizePlot",cellSizePlot(obj))
plotPDF(outDirPlot,dataset_name,patientID,"_02_GenesSizePlot",genesSizePlot(obj))
mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-")
plotPDF(outDirPlot,dataset_name,patientID,"_03_MitocondrialPlot",mit[["plot"]])

logThis("Raw data plotted")
###############################################################################
#'Removes the cells with the sum of genes counted for each cell higher than a 
#'certain threshold
i = 1 #counter
cells_count <- as.integer(args[3])

if(cells_count >= 0)
{
  #'If we inserted in the args the threshold, we remove them directly
  cells_to_rem <- getCells(obj)[getCellsSize(obj) > cells_count]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_04_CellSizePlot_cut_", i), cellSizePlot(obj))
} else {
  #'Otherwise, we ask for the threshold at runtime.
  #'It requires the input to be > 0
  while(cells_count < 0)
  {
    cells_count = as.integer(readline(prompt="Please insert a valid value for the cells count threshold:"))
    
    if(cells_count > 0)
    {
      cells_to_rem <- getCells(obj)[getCellsSize(obj) > cells_count]
      obj <- dropGenesCells(obj, cells = cells_to_rem)
      plotPDF(outDirPlot,dataset_name,patientID,paste0("_04_CellSizePlot_cut_", i),cellSizePlot(obj))
      
      answer = readline(prompt="Do you want to cut further? [y/n] ")
      if(answer == 'y')
        cells_count = -1
      
      i = i+1
    }
  }
}


###############################################################################
#'Removes the cells with with a number of genes with at least one read 
#'higher than a certain threshold
gene_count <- as.integer(args[4])
i = 1


if(gene_count >= 0)
{
  #'If we inserted in the args the threshold, we remove them directly
  columnCounts <- diff(obj@raw@p)
  names(columnCounts) <- obj@raw@Dimnames[[2]]
  cellGeneNumber <- sort(columnCounts, decreasing = FALSE) 
  cells_to_rem <- names(cellGeneNumber)[cellGeneNumber > gene_count]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_05_GeneCountPlot_cut_", i),genesSizePlot(obj))
} else{
  #'Otherwise, we ask for the threshold at runtime.
  #'It requires the input to be > 0
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
#'Removes the cells with with the percentage of mitochondrial higher than a
#'certain threshold

mitocondrial_count <- as.numeric(args[5])
i = 1

if(mitocondrial_count >= 0)
{
  #'If we inserted in the args the threshold, we remove them directly
  mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Analysis again
  to_rem <- mit[["sizes"]][["mit.percentage"]] > mitocondrial_count
  cells_to_rem <- rownames(mit[["sizes"]])[to_rem]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-")
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_06_MitocondrialCount_cut_", i),mit[["plot"]])
} else {
  while(mitocondrial_count < 0)
  {
    #'Otherwise, we ask for the threshold at runtime.
    #'It requires the input to be > 0
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

#'Removes the mitochondrial genes altogether
logThis("Removing the mitochondrial genes altogether")
genes_to_rem = getGenes(obj)[grep('^MT-', getGenes(obj))] 
cells_to_rem = getCells(obj)[which(getCellsSize(obj) == 0)]
obj = dropGenesCells(obj, genes_to_rem, cells_to_rem)

logThis("################################################################################")
logThis("Running clean")
#'Calls COTAN's clean() method
obj <- clean(obj)

logThis("Running cleanPlots")
cleanPlots <- cleanPlots(obj)

logThis("Data plotting")
plotPDF(outDirPlot,dataset_name,patientID,("_07_PCACells"),cleanPlots$pcaCells)
###############################################################################
#'Proceeds to remove (if asked) the "B" cluster
remove_clusterB <- as.integer(args[6])

#'if args[5] < 0, it means that the user wants to choose that specific moment
#'whether to remove the cluster or not by looking at the plot of the "nu" values
if(remove_clusterB < 0)
{
  answer = readline(prompt="Basing on the plots shown, do you want to remove the B cluster? [y/n] ")
  if(answer == 'y')
    remove_clusterB = 1
  else
    remove_clusterB = 0
}

#'if we want to remove the "B" cluster:
if(remove_clusterB == 1)
{
  #'Proceeds to remove it
  cells_to_rem <- rownames(cleanPlots$pcaCellsData)[cleanPlots$pcaCellsData[["groups"]] == "B"]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  
  #'Calls clean() and cleanPlots
  obj <- clean(obj)
  cleanPlots <- cleanPlots(obj)
  
  #'Useful plot
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_08_PCACellsBRemoval"),cleanPlots$pcaCells)
}

#'More plots
plotPDF(outDirPlot,dataset_name,patientID,"_09_CleanPlotGenes",cleanPlots$genes)
plotPDF(outDirPlot,dataset_name,patientID,"_10_CleanPlotUDE",cleanPlots$UDE)
plotPDF(outDirPlot,dataset_name,patientID,"_11_CleanPlotNu",cleanPlots$nu)
###############################################################################
#'Removes the low-efficiency cells with a low "nu" score

#'Plots useful to visualize the "nu" values
UDEPlot = cleanPlots$UDE
nuPlot = cleanPlots$nu

nu_threshold <- as.numeric(args[7])
i = 1

if(nu_threshold > 0)
{
  #'If we inserted in the args the threshold, we remove them directly
  nuDf = data.frame("nu" = sort(getNu(obj)), "n" = seq_along(getNu(obj)))
  cells_to_rem = rownames(nuDf)[nuDf[["nu"]] < nu_threshold]
  obj <- dropGenesCells(obj, cells = cells_to_rem)
  
  #'clean
  obj <- clean(obj)
  cleanPlots <- cleanPlots(obj)
  
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_12_PCACells_cut_", i),cleanPlots$pcaCells)
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_13_CleanPlotGenes_cut_", i),cleanPlots$genes)
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_14_CleanPlotUDE_cut_", i),cleanPlots$UDE)
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_15_CleanPlotNu_cut_", i),cleanPlots$nu)
} else {
  while(nu_threshold <= 0)
  {
    nu_threshold = as.numeric(readline(prompt="Please insert a valid value for the nu threshold: "))
    
    if(nu_threshold > 0)
    {
      nuDf = data.frame("nu" = sort(getNu(obj)), "n" = seq_along(getNu(obj)))
      cells_to_rem = rownames(nuDf)[nuDf[["nu"]] < nu_threshold]
      obj <- dropGenesCells(obj, cells = cells_to_rem)
      
      #'Otherwise, we ask for the threshold at runtime.
      #'It requires the input to be > 0
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
#'We can finally start the COTAN analysis. First of all, the DispersionBisection
logThis("Estimating the dispersion bijection")
obj = estimateDispersionBisection(obj, cores = 4)

#'And the coes matrices
logThis("Calculating the coex matrices (this might take a while)")
obj <- calculateCoex(obj)
################################################################################
logThis("GDI calculation")
#'Not mandatory, but we calculated the GDI anyway
quant.p = calculateGDI(obj)

#'Plots the GDI
GDIPlot = GDIPlot(obj, cond = paste0(" for dataset ", dataset_name, ", patient ", patientID), genes = knownCellsClean)
plotPDF(outDirPlot,dataset_name,patientID,"_16_GDIPlot",GDIPlot, width = 14, height = 14)
################################################################################
#' Uniform Clustering
logThis("Uniform clustering")
fineClusters <- cellsUniformClustering(obj, GDIThreshold = 1.4, cores = 6, saveObj = TRUE, outDir = outDirClustering)
obj <- addClusterization(obj, clName = "FineClusters", clusters = fineClusters)

#'Coex data for the clusters
logThis("Calculating the fine clusters' coex data")
c(coexDF, pValueDF) %<-% DEAOnClusters(obj, clusters = fineClusters)
obj <- addClusterizationCoex(obj, clName = "FineClusters", coexDF = coexDF)

logThis("Merging the clusters")
c(mergedClusters, coexDF, pValueDF) %<-%
  mergeUniformCellsClusters(obj, GDIThreshold = 1.4, cores = 6, saveObj = TRUE, outDir = outDirClustering)
obj <- addClusterization(obj, clName = "MergedClusters", clusters = mergedClusters, coexDF = coexDF)
################################################################################
plotPDF(outDirPlot,dataset_name,patientID,"_22_FineClustersSummary",clustersSummaryPlot(obj, clName = "FineClusters")[["plot"]])

#'Warning: if COTAN merges two or more clusters, the numbers in this
#'plots might be unreliable
plotPDF(outDirPlot,dataset_name,patientID,"_23_MergedClustersSummary",clustersSummaryPlot(obj, clName = "MergedClusters")[["plot"]])

#'Just a little extra: the clusters' markers (they mostly match the markers of
#'the known cell types regardless of whether we use the markers = knownCells
#'argument or not, which is why we have stuck to this method)
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

#'Prints the markers found
logThis("The known cell markers found in the data are: ")
markers <- cl
names(markers) <- gene
knownMarkersList = sort(markers[isMarker == 1])
knownMarkersList

###############################################################################
logThis("Plotting the cell type heatmap")
kCuts <- as.numeric(args[8])

#'kCuts is used only for visual reasons in the markers' heatmap
if(kCuts >= 1)
{
  #'If we inserted kCuts in the args the threshold, we calculate the heatmap directly
  clustersMarkersHeatmapPlot <- clustersMarkersHeatmapPlotB(
    obj,
    kCuts = kCuts,
    groupMarkers = knownCells,
    clName = "MergedClusters")[["heatmapPlot"]]
  plotPDF(outDirPlot,dataset_name,patientID,"_26_clustersMarkersHeatmapPlot",clustersMarkersHeatmapPlot, width=20, height=20)
} else {
  #'Otherwise, we ask for kCuts directly here, after plotting the "clean" heatmap
  clustersMarkersHeatmapPlot <- clustersMarkersHeatmapPlotB(
    obj,
    kCuts = 1,
    groupMarkers = knownCells,
    clName = "MergedClusters")[["heatmapPlot"]]
  plotPDF(outDirPlot,dataset_name,patientID,"_26_clustersMarkersHeatmapPlot",clustersMarkersHeatmapPlot, width=20, height=20)
  
  #'We need a number of cuts higher than 1
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

#'Computes the same heatmap calculated in the last section, only that this time
#'we need it so we can place it on the left of the enriched genes' heatmap as 
#'per the professor's instructions.
clustersMarkersPlot2 <- clustersMarkersHeatmapPlotB(obj, kCuts = kCuts,
                                                    groupMarkers = knownCells,
                                                    clName = "MergedClusters",
                                                    row_dend_width = 3.0,
                                                    cellsize = .2,
                                                    font_size = 6L,
                                                    use_cell_fun = FALSE)[["heatmapPlot"]]

if(enrichment_cut >= 1)
{
  #'If we have inserted the number of clusters manually, we plot it straight away
  c(enrichmentHm, enrichmentHmUnclustered, scoreDF, enrichmentClusters) %<-% EnrichmentHeatmap(obj, 
                                                                                               row_km = kCuts,
                                                                                               column_km = enrichment_cut,
                                                                                               groupMarkers = groupMarkers,
                                                                                               clName = "MergedClusters",
                                                                                               cellsHeatmap = clustersMarkersPlot2)
  names(enrichmentClusters) = as.character(seq(1, length(enrichmentClusters), by=1))
  save_list_to_csv(enrichmentClusters, outDir, "EnrichmentGenesClusters.csv", header = "cluster_ID,genes")
  
  plotPDF(outDirPlot,dataset_name,patientID,"_27_enrichmentHm",enrichmentHm, width = 35, height = 10)
} else {
  #Otherwise, we ask for the number of clusters after plotting the "unclustered" heatmap
  c(enrichmentHm, enrichmentHmUnclustered, scoreDF, enrichmentClusters) %<-% EnrichmentHeatmap(obj, row_km = kCuts,  column_km = 1, groupMarkers = groupMarkers, clName = "MergedClusters", cellsHeatmap = clustersMarkersPlot2)
  plotPDF(outDirPlot,dataset_name,patientID,"_27_enrichmentHm",enrichmentHm, width = 35, height = 10)
  
  #'We need at least one cluster
  while(enrichment_cut < 1)
  {
    enrichment_cut = as.numeric(readline(prompt="Please insert a valid value for the number of cuts: "))
    
    if(enrichment_cut >= 1)
    {
      c(enrichmentHm, enrichmentHmUnclustered, scoreDF, enrichmentClusters) %<-% EnrichmentHeatmap(obj, row_km = kCuts,  column_km = enrichment_cut, groupMarkers = groupMarkers, clName = "MergedClusters", cellsHeatmap = clustersMarkersPlot2)
      names(enrichmentClusters) = as.character(seq(1, length(enrichmentClusters), by=1))
      save_list_to_csv(enrichmentClusters, outDir, "EnrichmentGenesClusters.csv", header = "cluster_ID,genes")
      
      plotPDF(outDirPlot,dataset_name,patientID,"_27_enrichmentHm",enrichmentHm, width = 35, height = 10)
      
      answer = readline(prompt="Do you want to cut differently? [y/n]")
      if(answer == 'y')
        enrichment_cut = -1
    }
  }
}

plotPDF(outDirPlot,dataset_name,patientID,"_28_enrichmentHmUnclustered",enrichmentHmUnclustered, width = 35, height = 5)

################################################################################
logThis("Saving various clusters alternatives")

#'We calculate the same thing done in the last section but for every number
#'of clusters ranging between "start" and "end". Then, we save the results of the 
#'clusterization as a .csv file, and we save the relative heatmaps' .pdf(s) as well. 
#'This is useful when we need to calculate the optimal number of clusters for each donor
set.seed(1) #very important!!!
clustersMarkersPlot2 <- clustersMarkersHeatmapPlotB(obj, kCuts = kCuts,
                                                    groupMarkers = knownCells,
                                                    clName = "MergedClusters",
                                                    row_dend_width = 3.0,
                                                    cellsize = .2,
                                                    font_size = 6L,
                                                    use_cell_fun = FALSE)[["heatmapPlot"]]
start = 10
end = 20
for(i in start:end)
{
  set.seed(1)
  c(enrichmentHm, b, c, curr_cluster) %<-% EnrichmentHeatmap(obj, row_km = kCuts,  column_km = i, groupMarkers = groupMarkers, clName = "MergedClusters", cellsHeatmap = clustersMarkersPlot2)
  plotPDF(outDirPlot,dataset_name,patientID,paste0("_27_enrichmentHm_", i),enrichmentHm, width = 35, height = 10, plot_at_screen = FALSE)
  
  names(curr_cluster) = as.character(seq(1, length(curr_cluster), by=1))
  save_list_to_csv(curr_cluster, outDirEnrichmentCsv, paste0(dataset_name,"_",patientID,"_","EnrichmentGenesClusters_", i, ".csv"), header = "cluster_ID,genes")
}
################################################################################
rm(clustersMarkersPlot2)
gc()
################################################################################
#'Final logging
setLoggingLevel(newLevel = 2L) #makes the logging a little more sophisticated
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
