options(parallelly.fork.enable = TRUE)
library(Matrix)
library(data.table)
library(ggplot2)
library(ggrepel)
library(zeallot)
#devtools::install_github("seriph78/COTAN", ref = "devel")
library(COTAN)

setwd("~/Scrivania/CHL_Project")

dataset_name = "CO_STR"
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

#N104689 7189
#H197396 6548
#I130084 1356
patientID = "I130084"

cells_to_remove <- getCells(obj)[!grepl(patientID, colnames(mat))]
#cells_to_remove <- getCells(obj)[!(grepl("I130064|N130064|H158108" , colnames(mat)))]
#cells_to_remove <- getCells(obj)[!(grepl("I130064|H158108" , colnames(mat)))]

obj <- dropGenesCells(obj, cells = cells_to_remove)

################################################################################
#plots
ECDPlot(obj, yCut = 300)
cellSizePlot(obj) ##somma del numero di geni contati per ogni cellula (per ogni colonna (cellula), somma i valori righe (gene count)) 
genesSizePlot(obj) ##numeri di geni > 0
mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Percentuale di geni = MT (mitocondriali)
mit[["plot"]]

cells_to_rem <- getCells(obj)[getCellsSize(obj) > 10000]
obj <- dropGenesCells(obj, cells = cells_to_rem)
cellSizePlot(obj) ##somma del numero di geni contati per ogni cellula (per ogni colonna (cellula), somma i valori righe (gene count)) 

columnCounts <- diff(obj@raw@p)
names(columnCounts) <- obj@raw@Dimnames[[2]]
cellGeneNumber <- sort(columnCounts, decreasing = FALSE) 
cells_to_rem <- names(cellGeneNumber)[cellGeneNumber > 3000]
obj <- dropGenesCells(obj, cells = cells_to_rem)
genesSizePlot(obj)

mit <- mitochondrialPercentagePlot(obj, genePrefix = "^MT-") ##Analysis again
to_rem <- mit[["sizes"]][["mit.percentage"]] > 12.5
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
t <- dataset_name
outDir <- tempdir()
saveRDS(obj, file = file.path(outDir, paste0(t, ".cotan.RDS")))

################################################################################
quant.p = calculateGDI(obj)
head(quant.p)

#ottenuti dal paper

knownCells = list(
  "Activated fibroblasts CCL19+ ADAMADEC1+" 	= c('CCL19', 'ADAMDEC1'),
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

clustersSummaryPlot(obj, clName = "FineClusters")
clustersSummaryPlot(obj, clName = "MergedClusters")

c(cl, gene, score, pVal, adjPval, DEA, isMarker) %<-% findClustersMarkers(
  obj,
  n = 10L,
  clusters = mergedClusters,
  markers = knownCells,
  coexDF = coexDF,
  pValueDF = pValueDF,
  deltaExp = NULL,
  method = "bonferroni"
)

markers <- cl
names(markers) <- gene
markersInKnownCells = sort(markers[isMarker == 1])

plot1 <- UMAPPlot(t(obj@raw), clusters = fineClusters, title = "Fine Cluster UMAP Plot")
plot(plot1)


treePlot = clustersTreePlot(obj, kCuts = 17, clName = "MergedClusters" )
plot(treePlot$dend)
###############################################################################
library(rlang)
library(dendextend)
library(idyr)
library(grid) 
library(ComplexHeatmap)
library(circlize)

clustersMarkersHeatmapPlotB <- function(objCOTAN, groupMarkers, clName = NULL, kCuts = 3, conditionsList = NULL) 
{
  clName <- getClusterizationName(objCOTAN, clName = clName)
  
  expressionCl <- clustersDeltaExpression(objCOTAN, clName = clName)
  scoreDF <- geneSetEnrichment(groupMarkers = groupMarkers, clustersCoex = expressionCl)
  
  scoreDFT <- t(scoreDF[, 1L:(ncol(scoreDF) - 2L)])
  
  dend <- clustersTreePlot(objCOTAN, kCuts = kCuts)[["dend"]]
  
  dend = set(dend = dend, "branches_lwd", 2)
  {
    numDigits <- floor(log10(nrow(scoreDFT))) + 1L
    rownames(scoreDFT) <- formatC(as.numeric(rownames(scoreDFT)),
                                  width = numDigits, flag = "0")
  }
  
  if (!is_empty(conditionsList)) {
    # a data frame coming from clustersSummaryPlot()
    cond1 <- conditionsList[[1L]]
    
    if (is.numeric(cond1[["Cluster"]])) {
      numDigits <- floor(log10(length(cond1[["Cluster"]]))) + 1L
      cond1[["Cluster"]] <-
        formatC(cond1[["Cluster"]], width = numDigits, flag = "0")
    }
    
    cond1 <- cond1[, c("Cluster", "cond1" ,"CellNumber")] %>%
      pivot_wider(names_from = cond1, values_from = CellNumber)
    
    cond1[, 2L:3L] <- cond1[, 2L:3L] / rowSums(cond1[, 2L:3L])
    cond1 <- as.data.frame(cond1)
    rownames(cond1) <- cond1[["Cluster"]]
    cond1 <- cond1[rownames(scoreDFT), ]
    cond1[is.na(cond1)] <- 0
    
    # TODO: here instead of F and M we need a flexible number of conditions
    cond1Col <- c("F" = "deeppink1", "M" = "darkturquoise")
    hb <- rowAnnotation(
      cond1 = anno_barplot(cond1[, 2L:3L],
                           width = unit(3.0, "cm"),
                           gp = gpar(fill = cond1Col, col = "black"),
                           align_to = "right",
                           labels_gp = gpar(fontsize = 12L)),
      annotation_name_rot = 0L)
  } else {
    hb <- NULL
  }
  
  clsInfo <- clustersSummaryPlot(objCOTAN, clName = clName)[["data"]]
  {
    numDigits <- floor(log10(nrow(clsInfo))) + 1L
    clsInfo[["Cluster"]] <- formatC(clsInfo[["Cluster"]], width = numDigits, flag = "0")
  }
  
  rownames(clsInfo) <- clsInfo[["Cluster"]]
  clsInfo <- clsInfo[rownames(scoreDFT), ]
  
  freq <- set_names(clsInfo[["CellNumber"]], rownames(clsInfo))
  print(freq)
  freq2 <- set_names(paste0(clsInfo[["CellPercentage"]], "%"), rownames(clsInfo))
  
  #anno_numeric = anno_numeric(freq,bg_gp = gpar(fill = "orange", col = "black"))
  #ha <- rowAnnotation(cell.number = anno_numeric, annotation_name_rot = 0)
  
  ha2 <- rowAnnotation(cell.number = anno_text(freq2, gp = gpar(fontsize = 10)),
                       annotation_name_rot = 0)
  
  if (!is_empty(conditionsList)) {
    lgdList <- list(
      Legend(labels = c("Female", "Male"), title = "cond1",
             legend_gp = gpar(fill = cond1Col))
    )
  } else {
    lgdList <- list()
  }
  
  colorFunc <- colorRamp2(c(0, 1), c("lightblue", "red"))
  cellFunc <- function(j, i, x, y, width, height, fill) {
    grid.text(formatC(scoreDFT[i, j], digits = 1L, format = "f"),
              x, y, gp = gpar(fontsize = 9L))
  }
  
  finalHeatmap <- Heatmap(scoreDFT, rect_gp = gpar(col = "white", lwd = 1L),
                          cluster_rows = dend,
                          cluster_columns = FALSE,
                          col = colorFunc,
                          width = unit(28.0, "cm"),
                          row_dend_width = unit(8.0, "cm"),
                          #height = unit(6.0, "cm"),
                          column_names_gp = gpar(fontsize = 11L),
                          row_names_gp = gpar(fontsize = 11L),
                          cell_fun = cellFunc,
                          
                          #right_annotation = c(ha, ha2),
                          right_annotation = c(ha2),
                          
                          left_annotation = c(hb)
                          #, column_title = paste0("Gene marker set expression in ", sample)
  )
  
  finalHeatmap <- draw(finalHeatmap, annotation_legend_list = lgdList)
  
  return(list("heatmapPlot" = finalHeatmap, "dataScore" = scoreDF))
}

clustersMarkersHeatmapPlot(
  obj,
  kCuts = 17L,
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

EnrichmentHeatmap <- function(objCOTAN, groupMarkers, clName = NULL, kCuts = 3, conditionsList = NULL) 
{
  clName <- getClusterizationName(objCOTAN, clName = clName)
  
  expressionCl <- clustersDeltaExpression(objCOTAN, clName = clName)
  scoreDF <- geneSetEnrichment(groupMarkers = groupMarkers, clustersCoex = expressionCl)
  
  scoreDFT <- t(scoreDF[, 1L:(ncol(scoreDF) - 2L)])
  
  dend <- clustersTreePlot(objCOTAN, kCuts = kCuts)[["dend"]]
  
  dend = set(dend = dend, "branches_lwd", 2)
  {
    numDigits <- floor(log10(nrow(scoreDFT))) + 1L
    rownames(scoreDFT) <- formatC(as.numeric(rownames(scoreDFT)),
                                  width = numDigits, flag = "0")
  }
  
  if (!is_empty(conditionsList)) {
    # a data frame coming from clustersSummaryPlot()
    cond1 <- conditionsList[[1L]]
    
    if (is.numeric(cond1[["Cluster"]])) {
      numDigits <- floor(log10(length(cond1[["Cluster"]]))) + 1L
      cond1[["Cluster"]] <-
        formatC(cond1[["Cluster"]], width = numDigits, flag = "0")
    }
    
    cond1 <- cond1[, c("Cluster", "cond1" ,"CellNumber")] %>%
      pivot_wider(names_from = cond1, values_from = CellNumber)
    
    cond1[, 2L:3L] <- cond1[, 2L:3L] / rowSums(cond1[, 2L:3L])
    cond1 <- as.data.frame(cond1)
    rownames(cond1) <- cond1[["Cluster"]]
    cond1 <- cond1[rownames(scoreDFT), ]
    cond1[is.na(cond1)] <- 0
    
    # TODO: here instead of F and M we need a flexible number of conditions
    cond1Col <- c("F" = "deeppink1", "M" = "darkturquoise")
    hb <- rowAnnotation(
      cond1 = anno_barplot(cond1[, 2L:3L],
                           width = unit(3.0, "cm"),
                           gp = gpar(fill = cond1Col, col = "black"),
                           align_to = "right",
                           labels_gp = gpar(fontsize = 12L)),
      annotation_name_rot = 0L)
  } else {
    hb <- NULL
  }
  
  clsInfo <- clustersSummaryPlot(objCOTAN, clName = clName)[["data"]]
  {
    numDigits <- floor(log10(nrow(clsInfo))) + 1L
    clsInfo[["Cluster"]] <- formatC(clsInfo[["Cluster"]], width = numDigits, flag = "0")
  }
  
  rownames(clsInfo) <- clsInfo[["Cluster"]]
  clsInfo <- clsInfo[rownames(scoreDFT), ]
  
  freq <- set_names(clsInfo[["CellNumber"]], rownames(clsInfo))
  print(freq)
  freq2 <- set_names(paste0(clsInfo[["CellPercentage"]], "%"), rownames(clsInfo))
  
  #anno_numeric = anno_numeric(freq,bg_gp = gpar(fill = "orange", col = "black"))
  #ha <- rowAnnotation(cell.number = anno_numeric, annotation_name_rot = 0)
  
  ha2 <- rowAnnotation(cell.number = anno_text(freq2, gp = gpar(fontsize = 10)),
                       annotation_name_rot = 0)
  
  if (!is_empty(conditionsList)) {
    lgdList <- list(
      Legend(labels = c("Female", "Male"), title = "cond1",
             legend_gp = gpar(fill = cond1Col))
    )
  } else {
    lgdList <- list()
  }
  
  colorFunc <- colorRamp2(c(0, 1), c("lightblue", "red"))
  cellFunc <- function(j, i, x, y, width, height, fill) {
    #grid.text(formatC(scoreDFT[i, j], digits = 6L, format = "f"), x, y, gp = gpar(fontsize = 1L))
  }
  
  finalHeatmap <- Heatmap(scoreDFT, rect_gp = gpar(col = "white", lwd = 1L),
                          cluster_rows = dend,
                          cluster_columns = TRUE,
                          col = colorFunc,
                          width = unit(16, "cm"),
                          row_dend_width = unit(4.0, "cm"),
                          height = unit(4.0, "cm"),
                          column_names_gp = gpar(fontsize = 4L),
                          row_names_gp = gpar(fontsize = 4L),
                          cell_fun = cellFunc,
                          
                          #right_annotation = c(ha, ha2),
                          right_annotation = c(ha2),
                          
                          left_annotation = c(hb)
                          #, column_title = paste0("Gene marker set expression in ", sample)
  )
  
  finalHeatmap <- draw(finalHeatmap, annotation_legend_list = lgdList)
  
  return(list("heatmapPlot" = finalHeatmap, "dataScore" = scoreDF))
}

EnrichmentHeatmap(
  obj,
  groupMarkers = groupMarkers,
  clName = "MergedClusters" 
)
