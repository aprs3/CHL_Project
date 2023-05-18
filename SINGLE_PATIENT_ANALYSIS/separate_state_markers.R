setwd("~/Scrivania/Single_Patient_Analisys")
source("utils.R")

dataset = "CO_IMM"
S <- list("N130084", "I130084")


#to_load <- list(N, I, H)
states_to_load <- list(S)


p <- paste0(getwd(), "/", dataset)

optimal_clusters <- list()

i <- 1
cells_data  <- read_csv_data(paste0(p, "/known_cells_genes.csv"))
cells_names <- names(cells_data)

csvs <- list()

#TI_IMM_state_markers_I
for(i in 1:length(states_to_load))
{
  patient_type <- substr(states_to_load[[i]][[1]],1,1)
  
  #loads one of the known state markers file from the folder
  filename <- paste0(p, "/", "CO_IMM_I130084_N130084_OptimalIntersection.csv")
  curr_csv <- read_csv_data(filename)
  
  to_add <- c('cell_name,new_genes,known_genes')
  
  #for each cell we now need to get the newly found genes (meaning, the ones not
  #included in the known_cells_genes.csv file, which can be found by intersection)
  #they can be obtained simply by subtracting the known_cells_genes from it 
  for(i in 1:length(cells_data))
  {
    #loads the current cell's known genes
    known_cell_markers <- cells_data[[i]]
    #print(known_cell_markers)
    
    #current genes
    file_cell_marker <- curr_csv[[i]]
    
    #finds the new_genes by set subtraction, and the remainder are the already known genes
    new_genes           <- paste0(setdiff(file_cell_marker, known_cell_markers), collapse=';')
    already_known_genes <- paste0(intersect(file_cell_marker, known_cell_markers), collapse=';')
    
    tmp <- paste0(list(cells_names[[i]], new_genes, already_known_genes), collapse=',')
    to_add <- c(to_add, tmp)
  }
  
  path <- paste0(p, "/",dataset, "_state_markers_splitted_", patient_type, ".csv")
  file.create(path)
  file_conn = file(path)
  writeLines(to_add, file_conn)
  close(file_conn)
}


