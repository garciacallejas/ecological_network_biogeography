### Comparative biogeography of ecological networks

[![DOI](https://zenodo.org/badge/898574847.svg)](https://doi.org/10.5281/zenodo.14277893)

This repository holds a compendium of ecological networks of several interaction types, gathered from previous studies and databases, and harmonised. The repository also holds the code necessary to reproduce the findings of the associated study. For more information see the publication: 

https://doi.org/10.1101/2024.12.04.626839 

The dataset comes from different sources, either open or compiled by the authors, stored in the "data" folder and joined together in the following files:

- network metadata, including the ID of each network, spatial coordinates, and other info, is stored in the file "results/network_collection.csv"
- all nodes of all networks are stored in the file "results/network_nodes_collection.csv"
- all links of all networks are stored in the file "results/network_link_collection.csv"
  
#### Workflow

1) By running the scripts with prefix 00_*, data from each source is harmonised and stored in the "data" folder. Previous to this, a harmonised set of references is generated in the script with prefix 000.
2) Then, the whole dataset is joined together and filtered in the script with prefix 01.
3) Adjacency matrices and null network realisations are obtained in the scripts with prefix 02 and 03, and stored in the folders "data/(null_)(uni-bi)partite_adjacency_matrices".
4) Structural metrics from each network, and potentially from the null realisations, are obtained in the script with prefix 04 and stored in the "results/metrics" folder and the different "***null_(uni-bi)partite_metrics" folders.
5) Given the metric values from the original and null networks, z-scores are obtained in the script with prefix 05 and stored in "results/metrics_z_scores.csv".
6) The results and figures of the manuscript are obtained in scripts with prefix 06 and 07.

These scripts need a few R packages to run, see the relevant scripts for them. There are further auxiliary files with specific functions or supplementary information that are also included. The full dataset, results and figures are already generated at the relevant folders.

