#'This script works exactly like analyze_single_patients.R. However, it does not
#'extract the genes which do appear alongside a particular cell type's markers.
#'Instead, it works with .csv files which do already have stored, for each cell type, the
#'list of genes that appears to behave just like that cell type's markers. This is 
#'useful when we want, for example, to analyze which pathways are enriched by
#'the genes extracted by intersecting several clusters together from patients
#'with the same condition (refer to section 5.1 of the report for more details).
#'It is also suited for analyzing the pathways expressed by the "condition independent
#'genes) described in the report (section 8).
#'
#'Since the files that should be analyzed by this script do already contain
#'intersections of sets, there is no need to save them in a .csv. Instead, they
#'are saved in a .txt file which contains, alongside the name of the enriched pathway,
#'some useful information about the genes participating in the enrichment.
#'
#'usage: just set the "dataset" variable to the right folder. It should contain
#'a file named "known_cells_genes.csv" listing all the known cell types.

#'set it to whatever you like. The working directory path should contain at least
#'the analyze_intersection.R script
setwd("~/Scrivania/CHL_enrichment")
source("utils.R")

library(clusterProfiler)
library(limma)
library(ReactomePA)

#'Name of the folder containing the dataset (it should be inside the working directory)
dataset = "TI_STR_130064"

path = paste0(getwd(), "/", dataset, "/")
destination_path = paste0(path, "destinationFiles/")

dir.create(file.path(destination_path))

#gets the various files in the "dataset" folder
files <- list.files(path = path, pattern = "\\.csv$")  #loads all the .csv files in the dataset folder
files <- setdiff(files, list("known_cells_genes.csv")) #removes known_cells_genes.csv since we don't want to analyze it.

#'loads the cell data
known_cells_genes <- names(read_csv_data(paste0(path, 'known_cells_genes.csv')))

#for each file, we want to extract the pathways enriched:
for (file in files)
{
  #starts recording the print() output and saves it into a file named
  #[csv filename]_results.txt
  sink(paste0(destination_path, file, "_results.txt"))
  
  #loads the file
  file_data <- read_csv_data(paste0(path, file))
  
  #searches for the cell type in the file data
  for(cell in known_cells_genes)
  {
    #stores the cluster
    to_save <- file_data[[cell]]
    
    #removes cell names from the cluster. One way to do so would be to
    #subtract from it all strings appearing in known_cells_genes, like it's
    #a set operation
    to_save <- setdiff(to_save, known_cells_genes)
    
    #It might seem like we are saving a file, but we are not. This is just done
    #for verbose purposes.
    list_filename <- paste0(file, "_", cell, ".txt")
    dest <- paste0(destination_path, list_filename)
    print("============================================================================================================================")
    print(paste0("Cell Type: ", cell, ", source filename: ", file))
    
    #Converts the gene names into the ENTRZ ID format
    alias2Symbol(to_save, species = "Hs", expand.symbols = FALSE)
    to_save <- tryCatch(
      {
        bitr(to_save, fromType="SYMBOL",  toType = "ENTREZID", OrgDb = "org.Hs.eg.db") [["ENTREZID"]]
      },
      error=function(e) {
        message("Warning: the list contained no valid items")
        return(list())
      }
    )
    
    #if there are any "actual" genes in that cluster (and not just other cell types
    #which got clusterized alongside the current one)
    if(length(to_save) > 1)
    { 
      #DAVID (which searches for reactome, wiki pathways, and kegg)
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
        to_print <- dplyr::select(kkResult, Description, GeneRatio, pvalue, geneID, Count)
        to_print <- subset(to_print, pvalue <= 0.05)
        if(length(to_print) > 0)
        {
          curr_cell_pathways <- to_print$Description
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
    }
  }
  
  #print record stopped
  sink()
}

#safety measure
sink()