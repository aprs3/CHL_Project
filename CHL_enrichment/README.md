# Pathway enrichment analysis
We attempted to perform pathway enrichment analysis on the clusters of enriched genes obtained
by our experiments, in order to see whether we could discover any new pathway associated with the disease.
The procedure was the following. After we clusterized the enriched genes for each patient as described in section
5.1.1, we proceeded to extract the pathways associated to the genes which behaved like the markers of that particular
cell (meaning the ones that got clusterized alongside the markers themselves). We have used the clusterProfiler
package to query the KEGG database. The code to query Wiki Pathways, DAVID, and Reactome is left commented
for future studies.

## Scripts and file description
The files specific to the pathway enrichment analysis experiments were:
* **analyze single patients.R**: given the path to a folder containing a set of .csv files, each one listing the clusters
of enriched genes calculated for a specific patient as obtained by venn.R, it proceeds to extract, for each cell
listed in the known cells genes.csv file, the pathways enriched by the genes which do behave like the markers
of such cell.
* **analyze intersections.R**: like analyze single patients.R, but it calculates the pathway enrichment analysis on
files that do already list the genes behaving like a certain cluster, for each cell type.
* **venn.R**: a customized version of venn.R used which skips the search of the optimal clusterization for each
patient (thus assuming that it had already been found).
