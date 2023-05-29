#'Given the name of a folder containing a set of .csv files containing, each one,
#'a clusterization of the enrichment target genes for a particular patient, it  
#'proceeds to extract, for each cluster where the markers of a cell type listed  
#'in known_cells_genes.csv appear, the pathways enriched by the genes in that 
#'cluster. Finally, it saves the results both as a .csv file and as a .txt file as well.
#'
#'usage: just set the "dataset" variable to the right dataset. It should contain
#'a file named "known_cells_genes.csv" listing all the known cell types.

#'set it to whatever you like. The working directory path should contain at least
#'the analyze_single_patients.R script
setwd("~/Scrivania/CHL_enrichment")
source("utils.R")

library(clusterProfiler)
library(limma)
library(ReactomePA)

#'Name of the folder containing the dataset (it should be inside the working directory)
dataset = "CO_STR"

#'gets the dataset folder's path
path = paste0(getwd(), "/", dataset, "/")

#creates inside the dataset folder a subfolder that will host all the produced files.
destination_path = paste0(path, "destinationFiles/")
dir.create(file.path(destination_path))

files <- list.files(path = path, pattern = "\\.csv$")  #loads all the .csv files in the dataset folder
files <- setdiff(files, list("known_cells_genes.csv")) #removes known_cells_genes.csv from "files" since we don't want to analyze it.

#'loads the cell data
known_cells_genes <- names(read_csv_data(paste0(path, 'known_cells_genes.csv')))

#for each file, we want to extract the pathways enriched:
for (file in files)
{
  #loads the file
  file_data <- read_csv_data(paste0(path, file))
  to_save <- list()
  
  #starts recording the print() output, and saves it into a file named
  #[csv filename]_results.txt
  sink(paste0(destination_path, file, "_results.txt"))
  
  #searches for the cluster where a particular cell type appears
  for(cell in known_cells_genes)
  {
    curr_cell_pathways <- list()
    
    #searches for the cell type in the file data
    for(cluster in file_data)
    {
      #if the cell is in the currently analyzed cluster
      if(length(cluster) > 0 && cell %in% cluster)
      {
        #saves the cluster
        genes <- cluster
        
        #removes all cell names from the cluster. One way to do so would be to
        #subtract from it all strings appearing in known_cells_genes, like it's
        #a set operation
        genes <- setdiff(genes, known_cells_genes)
        
        #It might seem like we are saving a file, but we are not. This is just done
        #for verbose purposes.
        list_filename <- paste0(file, "_", cell, ".txt")
        dest <- paste0(destination_path, list_filename)
        print("===============================================================")
        print(paste0("Destination: ", dest))
        #print(genes)
        
        #Converts the gene names into the ENTRZ ID format
        alias2Symbol(genes, species = "Hs", expand.symbols = FALSE)
        genes <- bitr(genes, fromType="SYMBOL",  toType = "ENTREZID", OrgDb = "org.Hs.eg.db") [["ENTREZID"]]
        
        #DAVID (which searches for reactome, wiki pathways, and kegg)
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
          #gets the resulting dataframe
          kkResult <- kkResult@result
          
          #extracts the relevant columns
          to_print <- dplyr::select(kkResult, Description, GeneRatio, pvalue, geneID, Count)
          
          #for some reason, sometimes it still included the pathways with p-values too high,
          #so we remove them once more
          to_print <- subset(to_print, pvalue <= 0.05)
          if(length(to_print) > 0)
          {
            #stores the pathways' names (will be saved in the .csv file)
            #and prints the whole pathway data in the txt file.
            curr_cell_pathways <- to_print$Description
            print(to_print)
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
  #stops the recording of the print() output
  sink() 
  
  #saves, for each cell type, the list of pathways found.
  #filename: [dataset]_[patient_id]_clusterized_pathways.csv
  names(to_save) <- known_cells_genes
  
  save_list_to_csv(to_save, destination_path, paste0(substr(file, 1, 14), "_clusterized_pathways.csv"),
                   items_separator = '@', column_separator = '$', header = "cell_name$genes")
}

#safety measure
sink()