# Single patient analysis
Two patients (100064 and 100084) donated both inflamed and non-inflamed portions of their intestines in all the
examined datasets but CO STR. We were asked, as another extra, to calculate for each patient in each dataset the
intersection of the enrichment clusters obtained respectively from the inflamed and not inflamed portions of their
intestines, in order to obtain what we called ”condition independent genes”, meaning genes which do appear to
behave like the markers of a particular cell regardless of the condition of the intestine, which might be of independent
interest. We have used the same number of clusters for the enriched genes as calculated in section 5.1.1. The files
with the results can be found on the Github repository, inside the folder ”SINGLE PATIENT ANALYSIS”. (they
are the files named ”[Dataset name] [patient IDs] OptimalIntersection.csv”).


## Scripts and file description
To select manually the optimal number of cluster we used a custom version of find_state_markers.R.
