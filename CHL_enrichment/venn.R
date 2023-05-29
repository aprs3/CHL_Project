#'This script is a customized version of venn.R as found in the main folder which
#'skips the calculation of the optimal clusters entirely and just loads the 
#'<dataset>_<patientID>__clusterized_pathways.csv files present in the target folder
#'in order to extract, for each cell type, the pathways enriched by the genes
#'which appeared to behave constantly like the markers of that particular cell type
#'under a specific condition before intersecting them. This experiment is analogous
#'to the one described in section 5.1 but it intersects lists of pathways instead of lists of genes
#'
#'Usage: just edit the dataset and the to_load variables with, respectively, the
#'dataset name and the targeted patients' IDs.

#'set it to whatever you like. The working directory path should contain at least
#'the analyze_intersection.R script
setwd("~/Scrivania/CHL_enrichment")
library(VennDiagram)
library(zeallot)
source("utils.R")

#'Given a list of lists, it returns the intersection of those lists, as well as the
#'size of the union. Useful to calculate the Jaccard similarity score
lists_intersection <- function(l)
{
  curr_intersection <- l[[1]]
  curr_union <- l[[1]]
  
  #we could start from i = 2 but in case l has only one list it causes an exception
  for(i in 1:length(l))
  {
    curr_intersection <- intersect(curr_intersection, l[[i]])
    curr_union <- union(curr_union, l[[i]])
  }
  union_count <- length(curr_union)
  
  return(list("curr_intersection" <- curr_intersection, "union_count" <- union_count))
}

#dataset name. The folder should be placed inside the working directory.
dataset = "TI_IMM"

#patient IDs (this is just an example)
to_load <- c("I104689","I130064","I182231","I139892")
# to_load <- c("N109389","N119540","N130064","N158891")
# to_load <- c("H101694","H152638","H158108","H180844")

# dataset = "CO_IMM"
# to_load <- c("I114902","I121881","I130084","I175041")
# to_load <- c("N104689","N124246","N128400","N154787")
# to_load <- c("H139073","H197396")
# 
# dataset = "CO_STR"
# to_load <- c("I130084")
# to_load <- c("N107306","N124246","N104152","N104689")
# to_load <- c("H197396")
#  
# dataset = "TI_STR"
# to_load <- c("I104689","I130064")
# to_load <- c("N130064","N166301")
# to_load <- c("H158108","H180844")

#Get type of patient character for the folders
firstCharacter = substr(to_load[[1]],1,1)

#'Gets the character representing the condition of the patients 
#'('I' in the example above)
p <- paste0(getwd(), "/", dataset, "/destinationFiles")

#list containing the patients' data
patients <- list()

#loads the patients' data
for(patient in to_load)
{
  patient_filename <- paste0(p, "/", dataset, "_", patient, "_clusterized_pathways.csv")
  patients[[length(patients) + 1]] <- read_csv_data(patient_filename, csv_sep = '$', genes_sep = '@')
}

#loads the known cells' markers
cells_to_search <- names(read_csv_data(paste0(getwd(), "/", dataset, "/known_cells_genes.csv")))

intersections_list <- list()

#for each cell type we are targeting, we need to intersect the set of pathways
#found to be enriched by the genes which do behave like that cell's markers for
#the patients listed in to_load
for (cell_to_search in cells_to_search)
{
  print("#####################################################################")
  print(paste0("printing intersection for the cell ", cell_to_search))
  to_process <- list()
  for (patient in patients)
    to_process[[length(to_process) + 1]] <- patient[[cell_to_search]]
  
  #we calculate the intersection of such clusters
  c(combination_intersect, union_count) %<-% lists_intersection(to_process)
  intersections_list[[length(intersections_list) + 1]] <- combination_intersect
  
  print(combination_intersect)
}

#saves the data.
names(intersections_list) <- cells_to_search
filename = paste0(dataset, "_", paste(to_load, collapse = "_"), "_PathwaysIntersection.csv")
save_list_to_csv(intersections_list, paste0(p, "/"), filename,
                 items_separator = '@', column_separator = '$', header = "cell_name$genes")
