setwd("~/Scrivania/CHL_enrichment")
library(VennDiagram)
library(zeallot)
source("utils.R")

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

dataset = "TI_IMM"
to_load <- c("I104689","I130064","I182231","I139892")
to_load <- c("N109389","N119540","N130064","N158891")
to_load <- c("H101694","H152638","H158108","H180844")

dataset = "CO_IMM"
to_load <- c("I114902","I121881","I130084","I175041")
to_load <- c("N104689","N124246","N128400","N154787")
to_load <- c("H139073","H197396")

dataset = "CO_STR"
to_load <- c("I130084")
to_load <- c("N107306","N124246","N104152","N104689")
to_load <- c("H197396")
 
dataset = "TI_STR"
to_load <- c("I104689","I130064")
to_load <- c("N130064","N166301")
to_load <- c("H158108","H180844")

#Get type of patient character for the folders
firstCharacter = substr(to_load[[1]],1,1)

p <- paste0(getwd(), "/", dataset, "/destinationFiles")

#list containing the patients data
patients <- list()

for(patient in to_load)
{
  patient_filename <- paste0(p, "/", dataset, "_", patient, "_clusterized_pathways.csv")
  patients[[length(patients) + 1]] <- read_csv_data(patient_filename, csv_sep = '$', genes_sep = '@')
}

cells_to_search <- names(read_csv_data(paste0(getwd(), "/", dataset, "/known_cells_genes.csv")))

intersections_list <- list()
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

names(intersections_list) <- cells_to_search
filename = paste0(dataset, "_", paste(to_load, collapse = "_"), "_PathwaysIntersection.csv")
save_list_to_csv(intersections_list, paste0(p, "/"), filename,
                 items_separator = '@', column_separator = '$', header = "cell_name$genes")
