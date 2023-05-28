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

known_cells_genes <- names(read_csv_data(paste0(path, 'known_cells_genes.csv')))

for (file in files)
{
  #loads the file
  file_data <- read_csv_data(paste0(path, file))
  
  #searches for the cluster where a particular cell type resides
  for(cell in known_cells_genes)
  {
    to_save <- list()
    
    #searches for the cell type in the file data
    for(cluster in file_data)
    {
      #if the cell is in that cluster
      if(length(cluster) > 0 && cell %in% cluster)
      {
        #saves the cluster
        to_save <- cluster
        
        #removes cell names from the cluster. One way to do so would be to
        #subtract from it all strings appearing in known_cells_genes, like it's
        #a set operation
        to_save <- setdiff(to_save, known_cells_genes)
        
        list_filename <- paste0(file, "_", cell, ".txt")
        dest <- paste0(destination_path, list_filename)
        print("===============================================================")
        print(paste0("Destination: ", dest))
        #print(to_save)
        
        sink(dest)
        
        alias2Symbol(to_save, species = "Hs", expand.symbols = FALSE)
        to_save <- bitr(to_save, fromType="SYMBOL",  toType = "ENTREZID", OrgDb = "org.Hs.eg.db") [["ENTREZID"]]
        
        #DAVID (which searches for reactome, wiki pathways and kegg)
        #pathways <- list("EC_NUMBER" ,"REACTOME_PATHWAY" ,"WIKIPATHWAYS", "BBID", "BIOCARTA", "KEGG_PATHWAY")
        # pathways <- list("REACTOME_PATHWAY" ,"WIKIPATHWAYS","KEGG_PATHWAY")
        # for(pathway in pathways)
        # {
        #   print(pathway)
        #   david = enrichDAVID(gene = to_save, idType = "ENTREZ_GENE_ID",
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
        kkResult <- enrichKEGG(gene = to_save, organism = 'hsa', pvalueCutoff = 0.05)
        if(!is.null(kkResult))
        {
          kkResult <- kkResult@result
          to_print <- dplyr::select(kkResult, Description, GeneRatio, pvalue, geneID)
          to_print <- subset(to_print, pvalue <= 0.05)
          if(length(to_print) > 0)
          {
            print(to_print)
          }
        }
        
        #WIKI PATHWAY
        # enrichWPResult <- enrichWP(to_save, organism = "Homo sapiens", pvalueCutoff = 0.05)
        # if(!is.null(enrichWPResult))
        # {
        #   enrichWPResult <- enrichWPResult@result
        #   to_print <- dplyr::select(enrichWPResult, Description, GeneRatio, pvalue, geneID)
        #   to_print <- subset(to_print, pvalue <= 0.05)
        #   if(length(to_print) > 0)
        #   {
        #     print(to_print)
        #   }
        # }
        
        
        #REACTOME
        # enrichPathwayResult <- enrichPathway(gene=to_save, pvalueCutoff = 0.05, readable=TRUE)
        # if(!is.null(enrichPathwayResult))
        # {
        #   enrichPathwayResult <- enrichPathwayResult@result
        #   to_print <- dplyr::select(enrichPathwayResult, Description, GeneRatio, pvalue, geneID)
        #   to_print <- subset(to_print, pvalue <= 0.05)
        #   if(length(to_print) > 0)
        #   {
        #     print(to_print)
        #   }
        # }
        
        #lapply(to_save, write, dest, append=TRUE, ncolumns=1000)
        sink()
      }
    }
  }
}
