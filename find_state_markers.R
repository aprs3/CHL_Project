setwd("~/Scrivania/CHL_Project")
source("utils.R")
dataset = "TI_IMM"

N <- list("N109389", "N119540", "N130064", "N158891")
I <- list("I104689", "I130064", "I182231", "I139892")
H <- list("H101694", "H152638", "H158108", "H180844")

#to_load <- list(N, I, H)
states_to_load <- list(N, I, H)


p <- paste0(getwd(), "/", dataset)

optimal_clusters <- list()

i <- 1
cells_to_search <- names(read_csv_data(paste0(p, "/known_cells_genes.csv")))

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
  
  #for j in cell types
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

#TI_IMM_state_markers_I
for(i in 1:length(states_to_load))
{
  patient_type <- substr(states_to_load[[i]][[1]],1,1)
  filename <- paste0(dataset, "_state_markers_",patient_type, ".csv")
  
  save_list_to_csv(to_return[[i]], p, filename)
}
