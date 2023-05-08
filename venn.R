source("utils.R")
library(VennDiagram)

dataset = "TI_EPI"
to_load <- c("I115208", "I130064")

p <- paste0(getwd(), "/", dataset_name)

patients <- list()

for(patient in to_load)
{
  tmp <- read_csv_data(paste0(p, "/", patient, "/EnrichmentGenesClusters.csv"))
  patients[[length(patients)+1]] <- tmp
}

cell_to_search <- "Enterochromaffin cells"

to_process <- list()
for (patient in patients)
  for (cluster in patient)
    if(cell_to_search %in% cluster)
      to_process[[length(to_process) + 1]] <- cluster

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
filename <- paste0(p, "/", paste(to_load, collapse = "_"), "_venn.pdf")
venn.diagram(
  x = to_process,
  category.names = to_load,
  filename = filename
)

overlap <- calculate.overlap(x = to_process)
