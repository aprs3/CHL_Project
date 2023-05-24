# A Single Cell data analysis for Chron’s Disease with COTAN
2023 Computational Health Laboratory project.

## Group Members
- [Andrea Piras](https://github.com/aprs3)
- [Matteo Ninniri](https://github.com/Asduffo)

## Project description
This project makes single-cell RNA sequencing analysis with COTAN: In details it performs data cleaning, clustering of the cells, gene set enrichment. Given multiple patients analyzed can output for each cell the best intersection of genes between the patients. Last, finds the
exclusive genes for each condition (healthy, inflamed, not-inflamed samples) dividing them by new genes found in the enrichment and known genes from a given reference. 

## How to run
Clone the repository
```
git clone https://github.com/aprs3/CHL_Project.git
```
Install COTAN and his dependencies in R
```
devtools::install_github("seriph78/COTAN")
```

Download the [dataset](https://singlecell.broadinstitute.org/single_cell/study/SCP1884/human-cd-atlas-study-between-colon-and-terminal-ileum#study-download) files in the dataset folder. Compress with gzip the barcodes.tsv, features.tsv, raw.mtx. Then rename them features.tsv.zip, barcodes.tsv.zip, matrix.mtx.zip.

## Resources
- [COTAN](https://github.com/seriph78/COTAN/tree/devel)
- [Human CD atlas study between colon and terminal ileum](https://singlecell.broadinstitute.org/single_cell/study/SCP1884/human-cd-atlas-study-between-colon-and-terminal-ileum#study-summary)

