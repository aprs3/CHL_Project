library(rlang)
library(dendextend)
#library(idyr)
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
  
  finalHeatmap <- Heatmap(scoreDFT, 
                          rect_gp = gpar(col = "white", lwd = 1L),
                          cluster_rows = dend,
                          cluster_columns = TRUE,
                          col = colorFunc,
                          
                          row_dend_width = unit(3.0, "cm"),
                          
                          #width = unit(4.0, "cm"),
                          #height = unit(4.0, "cm"),
                          width = ncol(scoreDFT)*unit(.2, "cm"), 
                          height = nrow(scoreDFT)*unit(.2, "cm"),
                          
                          column_names_gp = gpar(fontsize = 6L),
                          row_names_gp = gpar(fontsize = 4L),
                          cell_fun = cellFunc,
                          
                          #right_annotation = c(ha, ha2),
                          #right_annotation = c(ha2),
                          
                          left_annotation = c(hb)
                          #, column_title = paste0("Gene marker set expression in ", sample)
  )
  
  finalHeatmapUnclustered <- Heatmap(scoreDFT, 
                          rect_gp = gpar(col = "white", lwd = 1L),
                          cluster_rows = dend,
                          cluster_columns = FALSE,
                          col = colorFunc,
                          
                          row_dend_width = unit(3.0, "cm"),
                          
                          #width = unit(4.0, "cm"),
                          #height = unit(4.0, "cm"),
                          width = ncol(scoreDFT)*unit(.2, "cm"), 
                          height = nrow(scoreDFT)*unit(.2, "cm"),
                          
                          column_names_gp = gpar(fontsize = 6L),
                          row_names_gp = gpar(fontsize = 4L),
                          cell_fun = cellFunc,
                          
                          #right_annotation = c(ha, ha2),
                          #right_annotation = c(ha2),
                          
                          left_annotation = c(hb)
                          #, column_title = paste0("Gene marker set expression in ", sample)
  )
  
  finalHeatmap <- draw(finalHeatmap, annotation_legend_list = lgdList)
  finalHeatmapUnclustered <- draw(finalHeatmapUnclustered, annotation_legend_list = lgdList)
  
  return(list("heatmapPlot" = finalHeatmap, 
              "heatmapPlotUnclustered" = finalHeatmapUnclustered, 
              "dataScore" = scoreDF))
} 
