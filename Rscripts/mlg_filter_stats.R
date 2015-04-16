#' # MLG Filter Statistics
#' 
#' ## Functions
#' 
#' See the source for details
#+ echo = FALSE
source("my_functions.R")
#'
#'
#' ### Directory for RAD seq data
RAD_data_dir <- "~/Documents/Grunwald/short-scripts/Genotype_error/"
#'
#' ### Analysis of *P. infestans*
library('poppr')
library('ape')
library('phangorn')
library('animation')
infdat <- RCurl::getURL("https://raw.githubusercontent.com/grunwaldlab/phytophthora_id/master/shiny-server/www/genoid_infestans/reduced_database.txt.csv")
infdat <- read.table(text = infdat, head = TRUE)
pinf <- df2genind(infdat[-c(1,2)], sep = "/", ploidy = 3, ind.names = infdat[[1]], pop = infdat[[2]])
ssr <- c(3,3,2,3,3,2,2,3,3,3,3,3)
x <- as.genclone(pinf)

fstats <- filter_stats(x, bruvo.dist, plot = TRUE, replen = ssr, loss = FALSE, nclone = 18)
title(main = expression(paste(italic("P. infestans"), " reference isolates (12 SSR loci)")))
# # 
# tiff(filename = "images/pinf_cluster.tiff", width = 85, height = 68, units = "mm", res = 1200, pointsize = 6)
# fstats <- filter_stats(x, bruvo.dist, plot = TRUE, replen = ssr, loss = FALSE, nclone = 18)
# title(main = expression(paste(italic("P. infestans"), " reference isolates (12 SSR loci)")))
# dev.off()

# legend("topright", legend = c("Nearest Neighbor", "UPGMA", "Farthest Neighbor"), 
#        col = c("red", "black", "blue"), pch = 1, title = "Clustering Method")
#' The plot above shows how multilocus genotypes collapse under differing 
#' algorithms over genetic distance. 
#' 
#' Below, we will collapse MLGs with a threshold of an average of 2 mutational
#' steps over all loci and create contingency tables relating the clustered
#' MLGs to the previously defined MLGs (eg. US-8).
z <- filter_stats(x, bruvo.dist, replen = ssr, loss = FALSE, threshold = 0.75/12, stats = "MLGS")
print(table(pop(x), z$farthest), zero.print = ".")
print(table(pop(x), z$nearest), zero.print = ".")
tab <- mlg.table(x, bar = FALSE)
colnames(tab) <- 1:ncol(tab)
print.table(tab, zero.print = ".") # contingency table for zero tolerance MLGs.
#'
#' Most of the MLGs were able to be resolved. US-8 and D.1 and D.2 are not
#' exactly resolved, but that is shown in the [trees produced](../images/upgma_average.html)
#' from the commented code below.
# Uncomment to regenerate plots
# 
# for (i in names(fstats)){
#   HTML_collapse(measure = i, x = x, treefun = "upgma", 
#                 distfun = "bruvo.dist", fstats = fstats, destdir = "images", 
#                 replen = ssr, loss = FALSE)  
#   HTML_collapse(measure = i, x = x, treefun = "nj", 
#                 distfun = "bruvo.dist", fstats = fstats, destdir = "images", 
#                 replen = ssr, loss = FALSE)  
# }

#' ### Simulated data 
#' 
#' This will create 20 populations with 20 samples and 10k SNPs. Each population
#' will have:
#' 
#'  - n = 20
#'  - 10,000 snps
#'  - an error rate of 0.1
#'  
#' In addition, half of these populations will have undergone one generation of
#' clonal reproduction. 
#' 
set.seed(20150415)
x <- lapply(1:10, getSims, n = 40, snps = 1e4, strucrat = 1, ploidy = 2, err = 0.1, na.perc = 0.21, clone = TRUE, n.cores = 4, mate_gen = 20)
y <- lapply(1:10, getSims, n = 40, snps = 1e4, strucrat = 1, ploidy = 2, err = 0.1, na.perc = 0.21, clone = FALSE, n.cores = 4, mate_gen = 20)
# x <- getSims(n = 200, snps = 1e4, strucrat = 1, ploidy = 2, err = 0.1, clone = TRUE, n.cores = 4)
# y <- getSims(n = 200, snps = 1e4, strucrat = 1, ploidy = 2, err = 0.1, clone = FALSE, n.cores = 4)
#'
#' For analysis, $\frac{1}{5}$th of the pooled samples will be kept.
z <- do.call("rbind", c(x, y))
z <- z[sample(nInd(z), nInd(z)/5)]
trueclones <- duplicated(substr(indNames(z), start = 1, stop = 10))
fstats <- filter_stats(z, bitwise.dist, plot = TRUE)
# the_threshold <- fstats$average$thresholds[sum(trueclones)] + .Machine$double.eps^0.5
the_threshold <- threshold_predictor(fstats$average$thresholds)
abline(v = the_threshold, lty = 2)
thresh     <- duplicated(mlg.filter(z, distance = bitwise.dist, 
                                    threshold = the_threshold, 
                                    algo = "a"))
(threshtable <- table(thresh, trueclones))
#' 
#' The tabulation is a power analysis of how many true and false positives there
#' are when collapsing at the threshold that gives the same number of known
#' clones/replicates.
#'
#' Below is labelling a tree with known clones.
#+ fig.width = 10, fig.height = 20
the_tree <- upgma(bitwise.dist(z))
clones <- substr(the_tree$tip.label[thresh], start = 1, stop = 10)
clones <- lapply(clones, grep, the_tree$tip.label)
edgelist <- length(which.edge(the_tree, the_tree$tip.label))
edgecols <- rep("black", edgelist)
for (i in clones){
  edgecols[which.edge(the_tree, the_tree$tip.label[i])] <- "red"
}
plot.phylo(the_tree, edge.color = edgecols, adj = 0, label.offset = 0.001)
axisPhylo(1)
title("Random sequences with 1000 SNPs and a 0.21 error rate")
#'
#+ 100reps, cache = TRUE, fig.show = "animate"
nreps <- 100
resarray <- array(data = integer(nreps*4), dim = c(2, 2, nreps), 
                  dimnames = c(dimnames(threshtable), NULL))
neararray <- resarray
fararray  <- resarray
avarray   <- resarray
samplist  <- lapply(1:nreps, function(x) list(samp = NULL, tree = NULL, 
                                              mlgs = NULL))
Sys.time()

for (i in 1:nreps){
  set.seed(i) # setting seed for accuracy.
  snps <- rpois(1, 1e3)
  samp1 <- lapply(1:10, getSims, n = 20, snps = snps, strucrat = 1, ploidy = 2, 
                  err = 0.05, na.perc = 0.21, clone = TRUE, n.cores = 4)
  samp2 <- lapply(1:10, getSims, n = 20, snps = snps, strucrat = 1, ploidy = 2, 
                  err = 0.05, na.perc = 0.21, clone = FALSE, n.cores = 4)
  samp <- do.call("rbind", c(samp1, samp2))
  samp@ploidy <- rep(2L, nInd(samp))
  samp <- samp[sample(nInd(samp), nInd(samp)/5)]
  trueclones <- duplicated(substr(indNames(samp), start = 1, stop = 10))
  fstats <- filter_stats(samp, bitwise.dist, plot = TRUE)
  # the_threshold <- fstats$average$thresholds[sum(trueclones)] + .Machine$double.eps^0.5
  title(paste("seed:", i, "n:", nInd(samp), "snps:", snps))
  the_threshold <- threshold_predictor(fstats$average$thresholds)
  the_distance  <- bitwise.dist(samp)
  z <- filter_stats(x = samp, distance = bitwise.dist, 
                    threshold = the_threshold, stats = "MLGs")
  abline(v = the_threshold, lty = 2)
  text(the_threshold, 0, 
       labels = paste("Threshold:", signif(the_threshold, 3)), 
       adj = 0)
  samplist[[i]]$samp <- samp
  samplist[[i]]$tree <- upgma(the_distance)
  samplist[[i]]$mlgs <- z
  
  athresh <- duplicated(z$average)
  nthresh <- duplicated(z$nearest)
  fthresh <- duplicated(z$farthest)
  
  avarray[, , i] <- table(athresh, trueclones)
  avarray[, , i] <- sweep(avarray[, , i], 2, colSums(avarray[, , i]), "/")

  neararray[, , i] <- table(nthresh, trueclones)
  neararray[, , i] <- sweep(neararray[, , i], 2, colSums(neararray[, , i]), "/")

  fararray[, , i] <- table(fthresh, trueclones)
  fararray[, , i] <- sweep(fararray[, , i], 2, colSums(fararray[, , i]), "/")
}
#' 
#' ### Results
#' 
Sys.time()
# color_mlg_tree(samp, upgma(bitwise.dist(samp)), z$average)
# axisPhylo(1)
#' 
#' Now we get to see how well we did.
(ares <- apply(avarray, 1:2, mean))
(nres <- apply(neararray, 1:2, mean))
(fres <- apply(fararray, 1:2, mean))
#' 
#+ results = 'asis', echo = FALSE
resmat <- matrix(signif(c(ares[2, 2], nres[2, 2], fres[2, 2], 
                        ares[2, 1], nres[2, 1], fres[2, 1])*100, 3),
                 ncol = 2, 
                 dimnames = list(method = names(z),
                                 c("True Positive %", "False Positive %")))
knitr::kable(resmat)
#' 
#' ### RAD seq data
#' 
#' Note that this data has no reference and has a lot of error. There are 10 
#' technical replicates. Each file represents a different parameter used for 
#' STACKS.
plinklist <- list(m3 = character(0), m4 = character(0), m10 = character(0), def = character(0))
plinklist[["m4"]]  <- "2R/PopSamples/data.out/PopSamples_m4/Popsouts_Rselec/out.replicates/plink.raw"
plinklist[["def"]] <- "2R/PopSamples/data.out/PopSamples_def/Popsouts_Rselec/out.replicates/plink.raw"
plinklist[["m3"]]  <- "2R/PopSamples/data.out/PopSamples_m3/Popsouts_Rselec/out.replicates/plink.raw"
plinklist[["m10"]] <- "2R/PopSamples/data.out/PopSamples_m10/Popsouts_Rselec/out.replicates/plink.raw"



contlist   <- plinklist # Contingency tables
threshlist <- plinklist # Threshold stats
#'
#' Steps:
#' 1. read in data
#' 2. mlg.filter on all algorithms and plot the thresholds.
#' 3. create the contingency table for each output (using UPGMA method).
#' 
#' Note for each plot regarding the MLG filter:
#' 
#' - Red: Nearest Neighbor clustering
#' - Blue: Farthest Neighbor clustering
#' - Black: UPGMA clustering (average neighbor)
#' 
#' The dotted lines represent the threshold at which the algorithms each
#' creates 10 clusters.
for (i in names(plinklist)){
  barb <- read.PLINK(paste(RAD_data_dir, plinklist[[i]], sep = "/"))
  show(barb)
  fstats <- filter_stats(barb, bitwise.dist, plot = TRUE, nclone = nInd(barb) - 10)
  title(paste(i, "SNPS:", nLoc(barb)))
  minthresh  <- fstats$average$thresholds[10] + .Machine$double.eps^0.5
#   nearthresh  <- fstats$nearest$thresholds[10] + .Machine$double.eps^0.5
#   farthresh  <- fstats$farthest$thresholds[10] + .Machine$double.eps^0.5
#   abline(v = minthresh, lty = 2)
#   abline(v = nearthresh, lty = 2, col = "red")
#   abline(v = farthresh, lty = 2, col = "blue")
#   legend("topright", legend = c("Nearest Neighbor", "UPGMA", "Farthest Neighbor"), 
#          col = c("red", "black", "blue"), pch = 1, title = "Clustering Method")
  thresh     <- mlg.filter(barb, distance = bitwise.dist, algorithm = "a", 
                           threshold = minthresh)
  trueclones <- vapply(strsplit(indNames(barb), "_"), "[[", character(1), 1)
  trueclones <- duplicated(trueclones)
  thresh     <- duplicated(thresh)
  contlist[[i]] <- table(thresh, trueclones)
  threshlist[[i]] <- fstats$average$thresholds
}
#'
#' Print the contingency tables and differences between threshold values to see
#' if there is a large jump indicating a separation between replicates and 
#' independent samples.
print(contlist)
for (i in threshlist){
  plot(diff(i), log = "y")
}

#' No such luck.
#' 
#' Now, we are plotting the tree where 10 samples are collapsed via average
#' neighbor (UPGMA) and color the tips with the true duplicates.
#' 
#' Note about the figure:
#' Tips are colored blue. Internal branches are colored black.
#' If the algorithm found samples with a distance below the threshold, their
#' connecting branches are colored red. The duplicated samples have red labels
#+ fig.width = 10, fig.height = 20
defupgma <- phangorn::upgma(bitwise.dist(barb))

z <- filter_stats(barb, bitwise.dist, threshold = minthresh, stats = "MLGS")
barbnames <- vapply(strsplit(indNames(barb), "_"), "[[", character(1), 1)
dupes <- barbnames[duplicated(barbnames)]
thecols <- ifelse(barbnames %in% dupes, "red", "black")
color_mlg_tree(barb, defupgma, z$average, tip.color = thecols)
axisPhylo(1)
#'
#' ## Session Info
#' 
options(width = 100)
devtools::session_info()