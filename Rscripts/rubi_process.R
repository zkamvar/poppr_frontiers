#' Quick and dirty analysis of *P. rubi* to detect technical reps 
#' ===================
#' 
#' First, we need to load the vcfR package and source the functions I have
#' written.
# devtools::install_github("knausb/vcfR@devel", build_vignettes = TRUE)
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
rubi_all <- create_chrom(name = "rubi", rubi)
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
library('poppr')
library('ape')
library('phangorn')
rubi.gt <- extract.gt(rubi_all, element = "GT", mask = TRUE)
rubi.gt[1:10, 1:10]
rubi.gt[rubi.gt == "./."] <- NA
rubi.gt[rubi.gt == "0/0"] <- 0
rubi.gt[rubi.gt == "1/1"] <- 2
rubi.gt[rubi.gt == "0/1"] <- 1
mode(rubi.gt) <- "integer"
rubi.gl <- new("genlight", t(rubi.gt), ind.names = colnames(rubi.gt))
samplenames <- vapply(strsplit(indNames(rubi.gl), "_"), "[", "a", 1) 
#'
#' Now, we can look at a tree:
#+ fig.width = 10, fig.height = 10
plot(rubi.gl)
#'
#+ fig.width = 5, fig.height = 10
plot(upgma(bitwise.dist(rubi.gl)))
axisPhylo(1)
#' 
#' The tree looks beautiful (kinda)!
#' 
#' Here's an example of filtering based on nearest neighbor distance.
#+ fig.width = 5, fig.height = 10
fs <- filter_stats(rubi.gl, distance = bitwise.dist, plot = TRUE)
(the_threshold <- threshold_predictor(fs$nearest$thresholds))
abline(v = the_threshold, lty = 2)
the_distance <- bitwise.dist(rubi.gl)
mlgs <- mlg.filter(rubi.gl, threshold = the_threshold, dist = the_distance, 
                   algorithm = "n")
color_mlg_tree(x = rubi.gl, tree = upgma(the_distance), newmlgs = mlgs)
axisPhylo(1)
#'
#' Admittedly, this looks... scary, but recalling that there is a LOT of missing
#' data, we might be able to trim out the isolates that have more than 25%
#' missing and get a better rate of success. 
allowed_missing <- round(nLoc(rubi.gl)*0.25)
miss <- vapply(rubi.gl@gen, function(x) length(x@NA.posi), integer(1))
rubi.nomiss <- rubi.gl[miss <= allowed_missing]
sample.nomiss <- samplenames[miss <= allowed_missing]
plot(rubi.nomiss)
#' 
#+ fig.width = 5, fig.height = 10
plot(upgma(bitwise.dist(rubi.nomiss)))
axisPhylo(1)
#' 
#' The tree looks beautiful (kinda)!
#' 
#' Here's an example of filtering based on nearest neighbor distance.
#+ fig.width = 5, fig.height = 10
fs <- filter_stats(rubi.nomiss, distance = bitwise.dist, plot = TRUE)
(the_threshold <- threshold_predictor(fs$average$thresholds, fraction = 0.75))
abline(v = the_threshold, lty = 2)
the_distance <- bitwise.dist(rubi.nomiss)
mlgs <- mlg.filter(rubi.nomiss, threshold = the_threshold, dist = the_distance, 
                   algorithm = "a")
color_mlg_tree(x = rubi.nomiss, tree = upgma(the_distance), newmlgs = mlgs)
axisPhylo(1)
table(duplicated(mlgs), duplicated(sample.nomiss))
#' 
#' ## Session Info
#' 
options(width = 100)
devtools::session_info()