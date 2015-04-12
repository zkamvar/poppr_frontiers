# R scripts for analysis

## mlg\_filter\_stats.R


Note that this script assumes the old version of adegenet. I have placed a new version of poppr adding in Jonah's changes with mine in this repository. 

This script's purpose is to explore how the MLG filtering methods perform in three ways:

 1. With *P. infestans* data from Genotype-ID
 2. With simulated data
 3. With RAD seq data from [DOI: 10.1111/1755-0998.12291](http://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12291/abstract)
 
### *P. infestans* setup. 
  packages: poppr, phangorn, ape, RCurl (you should have all of these), animate.

### Simulation setup

Same as previous

### RAD seq data setup

 1. Download the data from this [dryad link](http://datadryad.org/resource/doi:10.5061/dryad.g52m3). It is number 2.
 2. Change the path of RAD\_data\_dir in the R script.