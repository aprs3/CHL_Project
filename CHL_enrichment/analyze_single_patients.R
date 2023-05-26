setwd("~/Scrivania/CHL_enrichment")
source("utils.R")

library(clusterProfiler)
library(limma)
library(ReactomePA)

dataset = "TI_IMM"

path = paste0(getwd(), "/", dataset, "/")
destination_path = paste0(path, "destinationFiles/")

dir.create(file.path(destination_path))

#files <- list.files(path = path, pattern = to_load[[1]])
files <- list.files(path = path, pattern = "\\.csv$")
files <- setdiff(files, list("known_cells_genes.csv")) #removes known_cells_genes.csv since we don't want to analyze it.

known_cells_genes <- names(read_csv_data(paste0(path, 'known_cells_genes.csv')))

for (file in files)
{
  #loads the file
  file_data <- read_csv_data(paste0(path, file))
  to_save <- list()
  
  sink(paste0(destination_path, file, "_results.txt"))
  #searches for the cluster where a particular cell type resides
  for(cell in known_cells_genes)
  {
    curr_cell_pathways <- list()
    
    #searches for the cell type in the file data
    for(cluster in file_data)
    {
      #if the cell is in that cluster
      if(length(cluster) > 0 && cell %in% cluster)
      {
        #saves the cluster
        genes <- cluster
        
        #removes cell names from the cluster. One way to do so would be to
        #subtract from it all strings appearing in known_cells_genes, like it's
        #a set operation
        genes <- setdiff(genes, known_cells_genes)
        
        list_filename <- paste0(file, "_", cell, ".txt")
        dest <- paste0(destination_path, list_filename)
        print("===============================================================")
        print(paste0("Destination: ", dest))
        #print(genes)
        
        alias2Symbol(genes, species = "Hs", expand.symbols = FALSE)
        genes <- bitr(genes, fromType="SYMBOL",  toType = "ENTREZID", OrgDb = "org.Hs.eg.db") [["ENTREZID"]]
        
        #DAVID (which searches for reactome, wiki pathways and kegg)
        #pathways <- list("EC_NUMBER" ,"REACTOME_PATHWAY" ,"WIKIPATHWAYS", "BBID", "BIOCARTA", "KEGG_PATHWAY")
        # pathways <- list("REACTOME_PATHWAY" ,"WIKIPATHWAYS","KEGG_PATHWAY")
        # for(pathway in pathways)
        # {
        #   print(pathway)
        #   david = enrichDAVID(gene = genes, idType = "ENTREZ_GENE_ID",
        #                       annotation=pathway, david.user = "m.ninniri1@studenti.unipi.it",
        #                       pvalueCutoff = 0.05)@result
        #   to_print <- dplyr::select(david, Description, GeneRatio, pvalue, geneID)
        #   to_print <- subset(to_print, pvalue <= 0.05)
        #   if(length(to_print) > 0)
        #   {
        #     print(to_print)
        #   }
        # }
        
        #KEGG
        kkResult <- enrichKEGG(gene = genes, organism = 'hsa', pvalueCutoff = 0.05)
        if(!is.null(kkResult))
        {
          kkResult <- kkResult@result
          to_print <- dplyr::select(kkResult, Description, GeneRatio, pvalue, geneID)
          to_print <- subset(to_print, pvalue <= 0.05)
          if(length(to_print) > 0)
          {
            curr_cell_pathways <- to_print$Description
            print(curr_cell_pathways)
          }
        }
        
        #WIKI PATHWAY
        # enrichWPResult <- enrichWP(genes, organism = "Homo sapiens", pvalueCutoff = 0.05)
        # if(!is.null(enrichWPResult))
        # {
        #   enrichWPResult <- enrichWPResult@result
        #   to_print <- dplyr::select(enrichWPResult, Description, GeneRatio, pvalue, geneID)
        #   to_print <- subset(to_print, pvalue <= 0.05)
        #   if(length(to_print) > 0)
        #   {
        #     curr_cell_pathways <- to_print$Description
        #   }
        # }
        
        #REACTOME
        # enrichPathwayResult <- enrichPathway(gene=genes, pvalueCutoff = 0.05, readable=TRUE)
        # if(!is.null(enrichPathwayResult))
        # {
        #   enrichPathwayResult <- enrichPathwayResult@result
        #   to_print <- dplyr::select(enrichPathwayResult, Description, GeneRatio, pvalue, geneID)
        #   to_print <- subset(to_print, pvalue <= 0.05)
        #   if(length(to_print) > 0)
        #   {
        #     curr_cell_pathways <- to_print$Description
        #   }
        # }
        
        #lapply(genes, write, dest, append=TRUE, ncolumns=1000)
      }
    }
    
    to_save[[length(to_save) + 1]] <- as.list(curr_cell_pathways)
  }
  sink()
  
  names(to_save) <- known_cells_genes
  
  save_list_to_csv(to_save, destination_path, paste0(substr(file, 1, 14), "_clusterized_pathways.csv"),
                   items_separator = '@', column_separator = '$', header = "cell_name$genes")
}

sink()