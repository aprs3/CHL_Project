#'This script loads the .csv files containing, for each patient group, the intersection
#'for each cell type of the optimal clusters (as calculated by venn.R) containing
#'the genes which behaved like that specific cell's markers, and removes from each
#'patient group the union of the genes which appeared in the other two intersections.
#'The resulting files contain the genes which behave, for each cell type, like
#'that specific cell's markers only when the patient presents a specific condition
#'(inflamed, not inflamed, or healthy).
#'
#'Usage: just edit the "dataset" variable to the target dataset's name,
#'       and set the list of patient IDs "N", "I", and "H" with the respective IDs.

#'set it to whatever you like. The working directory path should contain, at least,
#'the main.R script
setwd("~/Scrivania/CHL_Project")
source("utils.R")

#Name of the folder containing the dataset (it should be inside the working directory)
dataset = "TI_IMM"

#list of patients' IDs for each condition
N <- list("N109389", "N119540", "N130064", "N158891")
I <- list("I104689", "I130064", "I182231", "I139892")
H <- list("H101694", "H152638", "H158108", "H180844")

states_to_load <- list(N, I, H)

#dataset's path
p <- paste0(getwd(), "/", dataset)


optimal_clusters <- list()

i <- 1
cells_to_search <- names(read_csv_data(paste0(p, "/known_cells_genes.csv")))

#loads the various <dataset name>_<patient IDs concatenation>__OptimalIntersection.csv files
for(state in states_to_load)
{
  filename = paste0(dataset, "_", paste(state, collapse = "_"), "_OptimalIntersection.csv")
  state_optimal_clusters <- read_csv_data(paste0(p, "/", filename), csv_sep = ',', genes_sep = ';')
  names(state_optimal_clusters) <- cells_to_search
  optimal_clusters[[length(optimal_clusters) + 1]] <- state_optimal_clusters
}

to_return <- list()

#for i in (number of patient types (H/I/N))
for(i in 1:length(optimal_clusters))
{
  to_return[[length(to_return) + 1]] <- list()
  
  #for each cell type
  for(j in 1:length(cells_to_search))
  {
    #quick init
    to_return[[i]][[length(to_return[[i]]) + 1]] <- list()
    
    #current cluster (we must subtract the union of the others patients from it)
    current_cluster <- optimal_clusters[[i]][[j]]
    
    remaining_patients_union <- list()
    
    #for each patient type (H/I/N)
    for(k in 1:length(optimal_clusters))
    {
      #'if this is not the current condition's intersection, it adds its elements
      #'to the current intersection 
      if(i != k)
      {
        remaining_patients_union <- union(remaining_patients_union, optimal_clusters[[k]][[j]])
      }
    }
    
    #performs the subtraction
    subtr <- setdiff(current_cluster, remaining_patients_union)
    
    to_return[[i]][[j]] <- subtr
  }
  
  names(to_return[[i]]) <- cells_to_search
}

#for each patient condition, saves the file with the new intersections as
#[dataset name]_state_markers_[patient condition].csv
for(i in 1:length(states_to_load))
{
  patient_type <- substr(states_to_load[[i]][[1]],1,1)
  filename <- paste0(dataset, "_state_markers_",patient_type, ".csv")
  
  save_list_to_csv(to_return[[i]], p, filename)
}
