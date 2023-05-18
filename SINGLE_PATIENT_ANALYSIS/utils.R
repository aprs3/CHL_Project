library(rlang)
library(dendextend)
#library(idyr)
library(grid) 
library(ComplexHeatmap)
library(circlize)

read_csv_data <- function(path, csv_sep = ',', genes_sep = ';')
{
  f = read.csv(path,sep=csv_sep)
  
  to_return <- c()
  
  for(genes in f$genes)
  {
    l <- strsplit(genes, genes_sep)
    to_return <- c(to_return, l)
  }
  
  names(to_return) <- f$cell_name
  return(to_return)
}

save_list_to_csv <- function(l, path, filename, header = "cell_name,genes")
{
  to_save <- c(header)
  i <- 1
  
  names_list <- names(l)
  
  for(e in l)
  {
    to_concat <- paste(l[[i]], collapse = ';')
    curr <- paste0(names(l)[[i]], ",", to_concat)
    
    to_save <- c(to_save, curr)
    
    i<-i+1
  }
  
  path <- paste0(path, "/", filename)
  file.create(path)
  file_conn = file(path)
  writeLines(to_save, file_conn)
  close(file_conn)
}

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

EnrichmentHeatmap <- function(objCOTAN, groupMarkers, clName = NULL, row_km = 1L, 
                              column_km = 1L, conditionsList = NULL,
                              row_dend_width = 3.0,
                              column_dend_height = 5.0,
                              cellsize = .2,
                              cellsHeatmap = NULL
                              ) 
{
  c(scoreDFT, dend, colorFunc, cellFunc,ha2, hb) %<-% clustersMarkersHeatmapPlot_main(objCOTAN, groupMarkers, clName, kCuts, conditionsList)
  
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
  
  clusters_data <- column_order(finalHeatmap)
  clusters_genes <- list()
  
  for (cluster in clusters_data) 
  {
    curr_cluster_genes <- c()
    
    for (cluster_member in cluster)
    {
      #print(names(groupMarkers)[[cluster_member]])
      curr_cluster_genes <- c(curr_cluster_genes, names(groupMarkers)[[cluster_member]])
    }
    
    print(curr_cluster_genes)
    curr_cluster_genes <- curr_cluster_genes[!duplicated(curr_cluster_genes)]
    
    to_add <- sort(curr_cluster_genes)
    clusters_genes[[length(clusters_genes)+1]] <- to_add
  }
  
  if(!is.null(cellsHeatmap))
    finalHeatmap <- cellsHeatmap + finalHeatmap
  
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
  
  #finalHeatmap <- draw(finalHeatmap, annotation_legend_list = lgdList)
  #finalHeatmapUnclustered <- draw(finalHeatmapUnclustered, annotation_legend_list = lgdList)
  
  return(list("heatmapPlot" = finalHeatmap, 
              "heatmapPlotUnclustered" = finalHeatmapUnclustered, 
              "dataScore" = scoreDFT,
              "clusters_genes" = clusters_genes))
} 


plotPDF <- function(path,dataset_name,patient_name,plotName, plt,plot_at_screen=TRUE, width = NULL, height = NULL) {
  if(plot_at_screen)
    plot(plt)
  
  pdf(file = paste0(path,dataset_name,"_",patient_name,"_",plotName,".pdf"),   # The directory you want to save the file in
      width = width, # The width of the plot in inches
      height = height) # The height of the plot in inches
  
  plot(plt)
  
  # Step 3: Run dev.off() to create the file!
  dev.off()
} 
