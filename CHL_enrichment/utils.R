#' This file contains several utilities to plot the clustering and gene enrichment data, 
#' saving the result of the enrichment as a .csv file, and reading such files.
#' It also features a slightly reworked clustersMarkersHeatmapPlot which partially 
#' fixes one bug in the COTAN library which might cause an exception if it 
#' merged two or more clusters in main.R. 

library(rlang)
library(dendextend)
library(grid) 
library(ComplexHeatmap)
library(circlize)

#'reads the data on a csv file and returns a named list where each entry is a list
#'itself containing the data of a particular row in the .csv. For example: if we do have
#'a file with the format
#'
#'cell,gene
#'CELL_1,GENE1;GENE2
#'CELL_2,GENE3;GENE4
#'
#'and we call read_csv_data(path, ',', ';')
#'it returns a list with the following format:
#'list({'CELL_1', list(GENE1, GENE2)},
#'     {'CELL_2', list(GENE3, GENE4)})
#'
#'csv_sep: the character which separates the row name (CELL_1 and CELL_2 in the
#'         example) from the elements in the list (anything on the left)
#'genes_sep: the character which splits the lists' elements
#'
#'The first row is always considered to be the header of the .csv file and is
#'ignored as a result. It is supposed to contain "cell,gene"
read_csv_data <- function(path, csv_sep = ',', genes_sep = ';')
{
  f = read.csv(path,sep=csv_sep)
  
  to_return <- c()
  
  #'For each row: 
  for(genes in f$genes)
  {
    #'creates the list of elements by splitting over genes_sep
    l <- strsplit(genes, genes_sep)
    to_return <- c(to_return, l)
  }
  
  #'names the elements in the list
  names(to_return) <- f$cell_name
  return(to_return)
}

#'saves a named list of lists as a .csv where each row is the content of one single
#'list (the first column is the name of the list and the second is the list of
#'items, split by items_separator. The two columns are separated by column_separator)
#'
#'example: if we have a named list of the type
#list({'CELL_1', list(GENE1, GENE2)},
#'    {'CELL_2', list(GENE3, GENE4)})
#'
#'calling save_list_to_csv(list, path, filename, "cell_name,genes", ';', ',')
#'saves the .csv file as:
#'cell,gene
#'CELL_1,GENE1;GENE2
#'CELL_2,GENE3;GENE4
save_list_to_csv <- function(l, path, filename, header = "cell_name,genes", items_separator = ';', column_separator = ",")
{
  to_save <- c(header)
  i <- 1
  
  names_list <- names(l)
  
  for(e in l)
  {
    to_concat <- paste(l[[i]], collapse = items_separator)
    curr <- paste0(names(l)[[i]], column_separator, to_concat)
    
    to_save <- c(to_save, curr)
    
    i<-i+1
  }
  
  path <- paste0(path, "/", filename)
  file.create(path)
  file_conn = file(path)
  writeLines(to_save, file_conn)
  close(file_conn)
}

#'This method is pretty much a copy of clustersMarkersHeatmapPlot() as seen in the
#'COTAN package. The only difference is that it does not compute the "ha" vector
#'which can sometimes cause some exceptions.
#'
#'parameters:
#'objCOTAN: the COTAN object
#'groupMarkers: list of markers
#'clName: name of the clusters to use (usually "MergedClusters")
#'kCuts: number of clusters for the heatmap rows
#'conditionsList: unused
clustersMarkersHeatmapPlot_main <- function(objCOTAN, groupMarkers, clName = NULL,
                                            kCuts = 3, conditionsList = NULL) 
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
  
  return(list("scoreDFT" <- scoreDFT, "dend" <- dend, "colorFunc" <- colorFunc,
              "cellFunc" <- cellFunc, "ha2" <- ha2, "hb" <- hb))
}

#'This method is similar to COTAN's clustersMarkersHeatmapPlot. The reason this 
#'portion of the code is not included in clustersMarkersHeatmapPlot_main defined
#'above is that the latter is also used by EnrichmentHeatmap defined below, 
#'so that this method only calls the final portion of clustersMarkersHeatmapPlot.
#'This makes the code more extensible.
#'
#'parameters:
#'objCOTAN: the COTAN object
#'groupMarkers: list of markers
#'clName = NULL: name of the clusters to be used (usually "MergedClusters")
#'kCuts = 3: number of clusters for the heatmap rows
#'conditionsList: unused
#'row_dend_width: dendrogram width
#'cellsize: the size of the cells squares
#'use_cell_fun: unused
#'font_size: font size
clustersMarkersHeatmapPlotB <- function(objCOTAN, 
                                        groupMarkers,
                                        clName = NULL,
                                        kCuts = 3,
                                        conditionsList = NULL,
                                        row_dend_width = 8.0, 
                                        cellsize = 1,
                                        use_cell_fun = TRUE,
                                        font_size = 11L)
{
  c(scoreDFT, dend, colorFunc,cellFunc,ha2, hb) %<-% clustersMarkersHeatmapPlot_main(objCOTAN,
                                                                                     groupMarkers,
                                                                                     clName,
                                                                                     kCuts,
                                                                                     conditionsList)
  right_annotation = c(ha2)
  left_annotation = c(hb)
  
  if(!use_cell_fun)
  {
    cellFunc = NULL
    right_annotation = NULL
    left_annotation= NULL
  }
  
  finalHeatmap <- Heatmap(scoreDFT, rect_gp = gpar(col = "white", lwd = 1L),
                          cluster_rows = dend,
                          cluster_columns = FALSE,
                          col = colorFunc,
                          row_dend_width = unit(row_dend_width, "cm"),
                          
                          width = ncol(scoreDFT)*unit(cellsize, "cm"), 
                          height = nrow(scoreDFT)*unit(cellsize, "cm"),
                          
                          column_names_gp = gpar(fontsize = font_size),
                          row_names_gp = gpar(fontsize = font_size),
                          cell_fun = cellFunc,
                          
                          right_annotation = right_annotation,
                          
                          left_annotation = left_annotation
  )
  
  #finalHeatmap <- draw(finalHeatmap, annotation_legend_list = lgdList)
  
  return(list("heatmapPlot" = finalHeatmap, "dataScore" = scoreDFT))
}

#'Calculates the heatmap which shows the enrichment of each gene passed as a 
#'parameter in the various clusters of the COTAN object passed in input.
#'It also returns the data of the various clusters so we can save them as a .csv file.
#'
#'objCOTAN: the COTAN object
#'groupMarkers: list of markers
#'clName: name of the clusters to be used (usually "MergedClusters")
#'row_km: number of row clusters
#'column_km: number of clusters to use for the target genes clusterization procedure
#'conditionsList: unused
#'row_dend_width: width of the row dendrogram
#'column_dend_height: height of the column dendrogram
#'cellsize: the size of each cell 
#'cellsHeatmap: a heatmap as calculated by clustersMarkersHeatmapPlotB. It will be
#'              placed on the left of the main heatmap.
#'
#'returns:
#'heatmapPlot": clustered heatmap plot
#'heatmapPlotUnclustered: unclustered heatmap plot
#'dataScore: raw data used to create the heatmap
#'clusters_genes: list of lists containing the various target enrichment genes' clusters
EnrichmentHeatmap <- function(objCOTAN, groupMarkers, clName = NULL, row_km = 1L, 
                              column_km = 1L, conditionsList = NULL,
                              row_dend_width = 3.0,
                              column_dend_height = 5.0,
                              cellsize = .2,
                              cellsHeatmap = NULL
                              ) 
{
  #important!
  set.seed(1)
  
  #Heatmap data
  c(scoreDFT, dend, colorFunc, cellFunc,ha2, hb) %<-% clustersMarkersHeatmapPlot_main(objCOTAN, groupMarkers, clName, kCuts, conditionsList)
  
  #the main heatmap
  finalHeatmap <- Heatmap(scoreDFT, 
                          rect_gp = gpar(col = "white", lwd = 1L),
                          cluster_rows = dend,
                          cluster_columns = TRUE,
                          col = colorFunc,
                          
                          row_dend_width = unit(row_dend_width, "cm"),
                          column_dend_height = unit(column_dend_height, "cm"),
                          
                          show_column_dend = TRUE,
                          column_km = column_km,
                          
                          width = ncol(scoreDFT)*unit(cellsize, "cm"), 
                          height = nrow(scoreDFT)*unit(cellsize, "cm"),
                          
                          column_names_gp = gpar(fontsize = 6L),
                          row_names_gp = gpar(fontsize = 4L),
                          
                          left_annotation = c(hb)
  )
  
  #Extracts the clusterization data
  clusters_data <- column_order(finalHeatmap)
  
  clusters_genes <- list()
  for (cluster in clusters_data) 
  {
    curr_cluster_genes <- list()
    
    for (cluster_member in cluster)
    {
      curr_cluster_genes[[length(curr_cluster_genes) + 1]] <- names(groupMarkers)[[cluster_member]]
    }
    
    to_add <- curr_cluster_genes
    
    clusters_genes[[length(clusters_genes)+1]] <- to_add
  }
  
  #merges the heatmap with cellsHeatmap, if we have passed it
  if(!is.null(cellsHeatmap))
    finalHeatmap <- cellsHeatmap + finalHeatmap
  
  #generates the unclustered heatmap
  finalHeatmapUnclustered <- Heatmap(scoreDFT, 
                          rect_gp = gpar(col = "white", lwd = 1L),
                          cluster_rows = dend,
                          cluster_columns = FALSE,
                          col = colorFunc,
                          
                          row_dend_width = unit(row_dend_width, "cm"),
                          
                          width = ncol(scoreDFT)*unit(cellsize, "cm"), 
                          height = nrow(scoreDFT)*unit(cellsize, "cm"),
                          
                          column_names_gp = gpar(fontsize = 6L),
                          row_names_gp = gpar(fontsize = 4L),
                          
                          left_annotation = c(hb)
  )
  
  
  if(!is.null(cellsHeatmap))
    finalHeatmapUnclustered <- cellsHeatmap + finalHeatmapUnclustered
  
  return(list("heatmapPlot" = finalHeatmap, 
              "heatmapPlotUnclustered" = finalHeatmapUnclustered, 
              "dataScore" = scoreDFT,
              "clusters_genes" = clusters_genes))
} 

#'plots the pdf passed and saves it as a .pdf file in the specified path.
#'
#'parameters:
#'path: working directory path
#'dataset_name: dataset name (es: TI_IMM)
#'patient_name: patient ID
#'plotName: file name
#'plt: the plot itself
#'plot_at_screen=TRUE: if true, display the plot as well on the screen
#'width = NULL: if not null, specifies the width of the image
#'height = NULL: if not null, specifies the height of the image
plotPDF <- function(path,dataset_name,patient_name,plotName, plt,plot_at_screen=TRUE, width = NULL, height = NULL) {
  if(plot_at_screen)
    plot(plt)
  
  pdf(file = paste0(path,dataset_name,"_",patient_name,"_",plotName,".pdf"),   # The directory you want to save the file in
      width = width,   # The width of the plot in inches
      height = height) # The height of the plot in inches
  
  set.seed(1) #very important since plot() re-calculates the heatmap
  plot(plt)
  
  dev.off()
} 
