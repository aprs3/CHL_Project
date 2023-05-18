setwd("~/Scrivania/Single_Patient_Analisys/")
library(VennDiagram)
library(zeallot)
source("utils.R")

lists_intersection <- function(l)
{
  curr_intersection <- l[[1]]
  curr_union <- l[[1]]
  
  for(i in 2: length(l))
  {
    curr_intersection <- intersect(curr_intersection, l[[i]])
    curr_union <- union(curr_union, l[[i]])
  }
  union_count <- length(curr_union)
  
  return(list("curr_intersection" <- curr_intersection, "union_count" <- union_count))
}

dataset = "CO_IMM"

to_load <- c("N130064", "I130064")

#to_load <- c("N104689","N154787","N128400","N124246")
#to_load <- c("H197396", "H139073")

#Get type of patient character for the folders
firstCharacter = substr(to_load[[1]],1,1)


p <- paste0(getwd(), "/", dataset)

#Creating clustering folder
outputDir <- paste0(getwd(), "/", dataset, "/venn_plots")
dir.create(file.path(outputDir))
outputDir <- paste0(getwd(), "/", dataset, "/venn_plots/", firstCharacter)
dir.create(file.path(outputDir))

start = 5
end = 20

#list containing for each patient the list of clusters (one for each cut)
#->[n_patiens, (start - end), dim_cluster]
patients <- list()

#list used to calculate all the combinations of cuts for all the patients
seqs <- list()

for(patient in to_load)
{
  #current patient list
  patients[[length(patients)+1]] <- list()
  seqs[[length(seqs)+1]]         <- seq.int(1, end - start + 1, 1)
  
  #puts into that list all the cuts data
  for(i in start:end)
  {
    tmp <- read_csv_data(paste0(p, "/", patient, "/","enrichment_csv/", dataset, "_", patient, "_", "EnrichmentGenesClusters_", i , ".csv"))
    patients[[length(patients)]][[i - (start - 1)]] <- tmp
  }
}

#number of combinations = (n cuts)^(n patients)
combinations <- expand.grid(seqs)

cells_to_search <- names(read_csv_data(paste0(p, "/known_cells_genes.csv")))

combination_average_score <- integer(length(combinations[[1]]))

#for each cell to search
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
      #print(paste0("cell_to_search = ", cell_to_search, ", curr_combination_idx = ", curr_combination_idx, ", curr_patient_idx = ", curr_patient_idx))
      
      #current patient
      curr_patient_cuts <- patients[[curr_patient_idx]]
      
      #current patient's cluster index
      curr_patient_cluster_idx <- combinations[[curr_patient_idx]][[curr_combination_idx]]
      
      #print(paste0("curr_patient_cuts ", curr_patient_cuts[[1]]))
      #print(paste0("curr_patient_cluster_idx ", curr_patient_cluster_idx[[140]]))
      
      #get the correct cluster data for the current patient
      curr_patient_correct_cluster <- curr_patient_cuts[[curr_patient_cluster_idx]]
      
      #print(paste0("curr_patient_correct_cluster ", curr_patient_correct_cluster))
      
      #extracts the cluster containing the target cell
      for (cluster in curr_patient_correct_cluster)
      {
        if(cell_to_search %in% cluster)
        {
          #print(paste0("cluster ", cluster))
          curr_combination[[length(curr_combination) + 1]] <- cluster
        }
      }
    }
    
    #now curr_combination is a list of length(to_load) elements where at position
    #i we have the cluster containing the target cell for the i^th patient.
    
    #we calculate the intersection of such clusters
    c(combination_intersect, union_count) %<-% lists_intersection(curr_combination)
    
    #calculates the score of the current combination, which is
    #(number of elements intersected)/(total number of elements in the clusters)
    curr_combination_score <- length(combination_intersect)/union_count
    
    combination_average_score[[curr_combination_idx]] <- combination_average_score[[curr_combination_idx]] + curr_combination_score
  }
}

combination_average_score <- combination_average_score/length(cells_to_search)

best_combination_index <- which.max(combination_average_score)
best_combination_value <- combination_average_score[[best_combination_index]]

best_combination_ncuts <- list()    #optimal number of clusters for each patient
best_combination_clusters <- list() #this contains the actual cuts for each patient 
for(i in 1:length(to_load))
{
  best_cut_for_patient_i <- combinations[[i]][[best_combination_index]] 
  best_combination_ncuts[[1+length(best_combination_ncuts)]] <- best_cut_for_patient_i + start - 1
  
  best_combination_clusters[[1+length(best_combination_clusters)]] <- patients[[i]][[best_cut_for_patient_i]]
}
best_combination_ncuts
best_combination_clusters

names(best_combination_ncuts) <- to_load
capture.output(best_combination_ncuts, file = paste0(p, "/", dataset, "_OptimalClusterCuts_", paste(to_load, collapse = "_") ,".txt"))
################################################################################
intersections_list <- list()
for (cell_to_search in cells_to_search)
{
  print("#####################################################################")
  print(paste0("printing intersection for the cell ", cell_to_search))
  to_process <- list()
  for (patient in best_combination_clusters)
    for (cluster in patient)
      if(cell_to_search %in% cluster)
        to_process[[length(to_process) + 1]] <- cluster
  
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
}

names(intersections_list) <- cells_to_search
filename = paste0(dataset, "_", paste(to_load, collapse = "_"), "_OptimalIntersection.csv")
save_list_to_csv(intersections_list, paste0(getwd(), "/", dataset, "/"), filename)
