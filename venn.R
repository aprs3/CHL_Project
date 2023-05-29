#' Given a list of patients, it calculates the optimal number of enriched genes
#' clusters for each patient (within a range set by the variables start and end)
#' such that the clusters where a cell type's markers are placed variates as
#' least as possible. Finally, it saves into a .txt file the optimal number of
#' clusters for each patient and the .csv file with the intersections found with
#' the method described in section 5.1.1. To execute it, set the variable "dataset" 
#' with the selected dataset name (for example TI_IMM), the list "to_load" with
#' the patients' IDs, and the variables start and end that specify the range 
#' where to search the optimal number of clusters.

#'set it to whatever you like. The working directory path should contain at least
#'the main.R script
setwd("~/Scrivania/CHL_Project")

#'various imports
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
to_load <- c("I104689", "I130064", "I182231", "I139892")

#'Gets the character representing the condition of the patients 
#'('I' in the example above)
firstCharacter = substr(to_load[[1]],1,1)

#'finds the dataset's path
p <- paste0(getwd(), "/", dataset)

#Creates the Venn plots folders inside the dataset's folder
outputDir <- paste0(getwd(), "/", dataset, "/venn_plots")
dir.create(file.path(outputDir))
outputDir <- paste0(getwd(), "/", dataset, "/venn_plots/", firstCharacter)
dir.create(file.path(outputDir))

#range start and end
start = 12
end = 20

#'A list containing for each patient the list of clusters (one for each clusters
#'number assignment)
#-> size: [n_patiens, (start - end), dim_cluster]
patients <- list()

#a list used to calculate all the combinations of cuts for all the patients
seqs <- list()

#'loads the patients' data
for(patient in to_load)
{
  #current patient list
  patients[[length(patients)+1]] <- list()
  seqs[[length(seqs)+1]]         <- seq.int(1, end - start + 1, 1)
  
  #puts into that list all the cuts' data
  for(i in start:end)
  {
    tmp <- read_csv_data(paste0(p, "/", patient, "/","enrichment_csv/", dataset, "_", patient, "_", "EnrichmentGenesClusters_", i , ".csv"))
    patients[[length(patients)]][[i - (start - 1)]] <- tmp
  }
}

#number of combinations = (n cuts)^(n patients)
combinations <- expand.grid(seqs)

#loads the known cells' markers
cells_to_search <- names(read_csv_data(paste0(p, "/known_cells_genes.csv")))

#this vector saves the average Jaccard score for each cluster number assignment among
#the various patients.
combination_average_score <- integer(length(combinations[[1]]))

#for each cell
for(cell_to_search in cells_to_search)
{
  print(cell_to_search)
  
  #for each combination of cluster cuts
  for (curr_combination_idx in 1:length(combinations[[1]]))
  {
    curr_combination <- list()
    
    #loads the combination's clusters
    for(curr_patient_idx in 1:length(to_load))
    {
      #current patient
      curr_patient_cuts <- patients[[curr_patient_idx]]
      
      #current patient cluster's index
      curr_patient_cluster_idx <- combinations[[curr_patient_idx]][[curr_combination_idx]]
      
      #get the correct cluster data for the current patient
      curr_patient_correct_cluster <- curr_patient_cuts[[curr_patient_cluster_idx]]
      
      #extracts the cluster containing the target cell
      for (cluster in curr_patient_correct_cluster)
      {
        if(cell_to_search %in% cluster)
        {
          curr_combination[[length(curr_combination) + 1]] <- cluster
        }
      }
    }
    
    #now curr_combination is a list of length(to_load) elements where at position
    #i we have the cluster containing the target cell for the i^th patient.
    
    #we calculate the intersection of such clusters
    c(combination_intersect, union_count) %<-% lists_intersection(curr_combination)
    
    #calculates the Jaccard of the current combination, which is
    #(number of elements intersected)/(total number of elements in the clusters)
    curr_combination_score <- length(combination_intersect)/union_count
    
    #'adds the score to the scores list
    combination_average_score[[curr_combination_idx]] <- combination_average_score[[curr_combination_idx]] + curr_combination_score
  }
}

#'calculates the average score for each assignment
combination_average_score <- combination_average_score/length(cells_to_search)

#'extracts the index of the optimal combination
best_combination_index <- which.max(combination_average_score)
best_combination_value <- combination_average_score[[best_combination_index]]

#'extracts the effective optimal cluster count's value for each patient
best_combination_ncuts <- list()    #optimal number of clusters for each patient
best_combination_clusters <- list() #this contains the actual cuts for each patient 


for(i in 1:length(to_load))
{
  best_cut_for_patient_i <- combinations[[i]][[best_combination_index]] 
  best_combination_ncuts[[1+length(best_combination_ncuts)]] <- best_cut_for_patient_i + start - 1
  
  best_combination_clusters[[1+length(best_combination_clusters)]] <- patients[[i]][[best_cut_for_patient_i]]
}

#quick print
best_combination_ncuts
best_combination_clusters

#names the lists and saves them into a file named [dataset]__OptimalClusterCuts_[list of patients in the vector to_load]
names(best_combination_ncuts) <- to_load
capture.output(best_combination_ncuts, file = paste0(p, "/", dataset, "_OptimalClusterCuts_", paste(to_load, collapse = "_") ,".txt"))
################################################################################

intersections_list <- list()

#saves the venn diagram for each cell type's optimal cluster for the patients in to_load
for (cell_to_search in cells_to_search)
{
  print("#####################################################################")
  print(paste0("printing intersection for the cell ", cell_to_search))
  to_process <- list()
  for (patient in best_combination_clusters)
    for (cluster in patient)
      if(cell_to_search %in% cluster)
        to_process[[length(to_process) + 1]] <- cluster
  
  #'now to_process contains the optimal clusters
  #'for a particular cell type (one cluster for each patient in to_load)
  
  #'plots the venn diagram
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  filename <- paste0(outputDir, "/", paste(to_load, collapse = "_"), "_venn_", cell_to_search, ".pdf")
  venn.diagram(
    x = to_process,
    category.names = to_load,
    print.mode=c("raw","percent"),
    filename = filename
  )
  
  #we calculate the intersection of such clusters
  c(combination_intersect, union_count) %<-% lists_intersection(to_process)
  intersections_list[[length(intersections_list) + 1]] <- combination_intersect
  
  print(combination_intersect)
}

#'saves the intersection of the various patients' optimal clusters (one intersection
#'for each cell type), in a file named
#'<dataset name>_<concatenations of the patients' IDs>_OptimalIntersection.csv
names(intersections_list) <- cells_to_search
filename = paste0(dataset, "_", paste(to_load, collapse = "_"), "_OptimalIntersection.csv")
save_list_to_csv(intersections_list, paste0(getwd(), "/", dataset, "/"), filename)
