#' Quick and dirty analysis of *P. rubi* to detect technical reps 
#' ===================
#' 
#' First, we need to load the vcfR package and source the functions I have
#' written.
devtools::install_github("knausb/vcfR@devel", build_vignettes = TRUE)
library("vcfR")
source("my_functions.R")
#' Now we read in the rubi data. Note that you must have the data in the folder
#' above for this to work.
rubi <- read.vcf("../All_rubi.vcf", limit = 1e+10)
#'
#' ## Cleaning
#' 
#' Since the rubi data is currently unfiltered, we will use Brian's methods
#' to filter it.
#' 
#' I attempted to utilize scaffold 1, but it had < 200 variants, so I'm
#' reading in the whole damn thing. 
#+ message = FALSE
rubi_all <- create_chrom(name = "rubi", rubi1)
#'
#' Now to plot the depth info. According to Brian's manual, this is necessary
#' to determine if the depth in the INFO column is the same as the depth reported
#' with the genotypes.
dp <- extract.gt(rubi_all, element = "DP", as.numeric = TRUE)
plot(rowSums(dp), rubi_all@var.info$DP, xlab = "Depth reported with genotypes", 
     ylab = "Depth reported in INFO field")
abline(a=0, b=1)
#'
#' Since they are not the same, we will change the depth to the depth
#' reported with the genotypes
rubi_all@var.info$DP <- rowSums(dp)
#'
plot(rubi_all)
#'
#' I'm kinda just shooting in the dark with this.
options(width = 100)
rubi_all <- masker(rubi_all, min_MQ = 43.9, max_MQ = 61, min_QUAL = 999, 
                   min_DP = 1000)
plot(rubi_all)
head(rubi_all)
#'
#' ## Transferring to genlight object
#' 
rubi.gt <- extract.gt(rubi_all, element = "GT")
rubi.gt[1:10, 1:10]
rubi.gt[rubi.gt == "./."] <- NA
rubi.gt[rubi.gt == "0/0"] <- 0
rubi.gt[rubi.gt == "1/1"] <- 2
rubi.gt[rubi.gt == "0/1"] <- 1
rubi.gt <- rubi.gt[rubi_all@var.info$mask, ]
mode(rubi.gt) <- "integer"
rubi.gl <- new("genlight", t(rubi.gt), ind.names = colnames(rubi.gt))
#'
#' Now, we can look at a tree:
library('poppr')
library('ape')
library('phangorn')
plot(upgma(bitwise.dist(rubi.gl)))
axisPhylo(3)
#' ## Session Info
#' 
options(width = 100)
devtools::session_info()