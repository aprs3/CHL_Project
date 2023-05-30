# Pathway enrichment analysis
We attempted to perform pathway enrichment analysis on the clusters of enriched genes obtained by our experiments, in order to see whether we could discover any new pathway associated with the disease. 

We have used the clusterProfiler package to query the KEGG database. The code to query Wiki Pathways, DAVID, and Reactome is left commented for future studies.

## Folders description
### CO_IMM, TI_IMM, CO_STR, and TI_STR
This folders contains the results of the pathway enrichment analysis on the optimal clusters found for each patient with a particular condition in each dataset.
Each folder contains a set of .csv files with a filename with the format <dataset name>_<patient ID>_EnrichmentGenesClusters_<number of clusters of genes>.csv, where "number of clusters of genes" is the optimal number of clusters found by the main experiment for that particular patient.

Finally, each folder contains a subfolder named "destinationFiles". Inside the folders there are three types of files: 
* **<dataset name>_<patient ID>_EnrichmentGenesClusters_<number of clusters of genes>.csv_results.txt**: It contains, for each cell type present in known_cells_genes.csv, the pathways enriched by the genes which appeared in the cluster where that cell type was clusterized inside the file <dataset name>_<patient ID>_EnrichmentGenesClusters_<number of clusters of genes>.csv.
* **<dataset name>_<patient ID>_clusterized_pathways.csv**:it contains the pathways listed in <dataset name>_<patient ID>_EnrichmentGenesClusters_<number of clusters of genes>.csv_results.txt, but listed as a .csv file. The cell column is separated by the pathways column with a "$", while the pathways themselves are separated by a "@".
* **<dataset name>_<list of patient IDs>_PathwaysIntersection.csv**: It contains, for each cell type in known_cells_genes.csv, the intersection of the pathways present in that cell in the files <dataset name>_<patient ID>_clusterized_pathways.csv, where the patient ID is in the list of patient IDs in the filename. The cell column is separated by the pathways column with a "$", while the pathways themselves are separated by a "@".

### CO_IMM_intersect, TI_IMM_intersect, CO_STR_intersect, and TI_STR_intersect
This folders contains the results of the pathway enrichment analysis performed on the intersections of the optimal clusters found for each patient with a particular condition in each dataset (this means the files present in the main experiments named <dataset name>_<list of patient IDs>_OptimalIntersection.csv, of which the relevant copies are placed in this folders). The "destinationFiles" folder contains only one type of files: <dataset name>_<list of patient IDs>_OptimalIntersection.csv_results.txt, which lists the pathways found for each cell type listed in <dataset name>_<list of patient IDs>_OptimalIntersection.csv, the pathways associated to the genes listed in the same row in the .csv file.

### CO_IMM_130084, TI_IMM_130064, and TI_STR_130064


## Scripts and file description
The files specific to the pathway enrichment analysis experiments were:
* **analyze single patients.R**: given the path to a folder containing a set of .csv files, each one listing the clusters of enriched genes calculated for a specific patient as obtained by venn.R, it proceeds to extract, for each cell listed in the known cells genes.csv file, the pathways enriched by the genes which do behave like the markers of such cell.
* **analyze intersections.R**: like analyze single patients.R, but it calculates the pathway enrichment analysis on files that do already list the genes behaving like a certain cluster, for each cell type.
* **venn.R**: a customized version of venn.R used which skips the search of the optimal clusterization for each patient (thus assuming that it had already been found).
