#' # MLG Filter Statistics
#' 
#' ## Functions
#' 
#' See the source for details
#+ echo = FALSE
getNames <- function(n){
  vapply(1:n, function(x) paste(sample(letters, 10), collapse = ""), character(1))
}

getSims <- function(z = 1, n = 10, snps = 1e6, strucrat = c(0.25, 0.75), 
                    clone = TRUE, err = 0.1, n.cores = 4, ...){
  lam <- sample(strucrat, 1)
  if (lam == 1){
    perc <- 1
  } else {
    perc <- ceiling(sample(rpois(1000, lambda = snps*lam), 1)/snps)    
  }

  n <- sample(rpois(1000, lambda = n), 1)
  the_names <- getNames(n)
  res <- glSim(n, n.snp.nonstruc = snps*perc, n.snp.struc = snps*(1-perc), n.cores = n.cores, ...)
  if (clone){
    clones <- sample(n, replace = TRUE)
    the_names <- paste(the_names[clones], 1:n, sep = ".")
    res <- res[clones]
    clones <- duplicated(clones)
  }
  mat   <- as.matrix(res)
  nerrs <- round(ncol(mat)*err)
  for (i in seq(nrow(mat))){
    mat[i, sample(ncol(mat), nerrs)] <- sample(c(0:2, NA), nerrs, replace = TRUE)
  }
  if (clone){
    for (i in clones){
      mat[i, sample(ncol(mat), nerrs)] <- sample(c(0:2, NA), nerrs, replace = TRUE)      
    }
  } 

  res <- new("genlight", mat, ploidy = 2, ind.names = the_names, n.cores = n.cores)
  return(res)
}


filter_stats <- function(x, distance = bitwise.dist, threshold = 1, 
                         stats = "All", missing = "ignore", plot = FALSE, nclone = NULL, ...){
  distmat <- distance(x, ...)
  f <- mlg.filter(x, threshold, missing, algorithm = "f", distance = distmat, 
                  stats = stats, ...)
  a <- mlg.filter(x, threshold, missing, algorithm = "a", distance = distmat, 
                  stats = stats, ...)
  n <- mlg.filter(x, threshold, missing, algorithm = "n", distance = distmat, 
                  stats = stats, ...)
  fanlist <- list(farthest = f, average = a, nearest = n)
  if (stats == "All"){
    for (i in names(fanlist)){
      names(fanlist[[i]]) <- c("MLG", "thresholds", "mat", "size")
    }
    if (plot){
      plot_filter_stats(x, fanlist, distmat, nclone)
    }
  }
  return(fanlist)
}

plot_filter_stats <- function(x, fstats, distmat, nclone = NULL){
  upper <- round(max(distmat), digits = 1)
  ylims <- c(ifelse(is.genind(x), mlg(x, quiet = TRUE), nInd(x)), 1)
  plot(x = c(upper, 0), y = ylims, type = "n",
       ylab = "Number of Multilocus Lineages",
       xlab = "Genetic Distance Cutoff")
  a <- fstats$average$thresholds
  n <- fstats$nearest$thresholds
  f <- fstats$farthest$thresholds
  points(x = rev(a), y = 1:length(a))
  points(x = rev(f), y = 1:length(f), col = "blue")
  points(x = rev(n), y = 1:length(n), col = "red") 
  if (!is.null(nclone)){
    abline(v = a[1 + length(a) - nclone] + .Machine$double.eps^0.5, lty = 2)
    abline(v = f[1 + length(f) - nclone] + .Machine$double.eps^0.5, lty = 2, 
           col = "blue")
    abline(v = n[1 + length(n) - nclone] + .Machine$double.eps^0.5, lty = 2, 
           col = "red")
    abline(h = nclone)
    text(upper, nclone, labels = paste0("n = ", nclone), 
         adj = c(1, -0.5))
    legend("topright", 
           legend = c("Nearest Neighbor", "UPGMA", "Farthest Neighbor"), 
           col = c("red", "black", "blue"), 
           pch = 1, 
           lty = 2,
           title = "Clustering Method")
  } else {
    legend("topright", 
           legend = c("Nearest Neighbor", "UPGMA", "Farthest Neighbor"), 
           col = c("red", "black", "blue"), 
           pch = 1, 
           title = "Clustering Method")    
  }
}

color_mlg_tree <- function(x, tree, newmlgs, ...){
  mlg_vector <- newmlgs
  mlg_cols <- "red"
  all_edges <- length(which.edge(tree, tree$tip.label))
  edge_cols <- rep("black", all_edges)
  for (i in tree$tip.label){
    edge_cols[which.edge(tree, i)] <- "blue"
  }
  edge_widths <- rep(1, all_edges)
  for (i in unique(mlg_vector)){
    indices <- which.edge(tree, tree$tip.label[mlg_vector == i])
    if (length(indices) > 1){
      edge_cols[indices] <- "red"      
      edge_widths[indices] <- 3
    }
  }
  tree$tip.label <- paste(tree$tip.label, pop(x), sep = "_")
  plot.phylo(tree, edge.color = edge_cols, edge.width = edge_widths, adj = 0, 
             label.offset = 0.005, ...)
}

HTML_collapse <- function(measure, x, treefun, distfun, fstats, destdir = NULL, ...){
  desc <- "An exploration of collapsing MLGs using"
  descname <- switch(measure, nearest = "Nearest Neighbor",
                     farthest = "Farthest Neighbor",
                     average = "UPGMA")
  desc <- paste(desc, descname)
  TREEFUN <- match.fun(treefun)
  DISTFUN <- match.fun(distfun)
  cwd <- getwd()
  if (!is.null(destdir)){
    setwd(destdir)
  }
  the_tree <- TREEFUN(DISTFUN(x, ...))
  isnj <- length(grep("nj", substitute(treefun))) > 0
  if (isnj){
    the_tree <- ladderize(poppr:::fix_negative_branch(the_tree))
  }
  saveHTML({
    for (i in unique(fstats[[measure]]$thresholds)){
      z <- mlg.filter(x, threshold = i, missing = "ignore", algorithm = measure, 
                      distance = distfun, stats = "MLGs", ...)
      if (isnj){
        color_mlg_tree(x, the_tree, z, type = "u", lab4ut = "axial")
      } else {
        color_mlg_tree(x, the_tree, z)
      }
      
      if (isnj){
        add.scale.bar()
      } else {
        axisPhylo(1)
      }
      title(paste("Threshold:", round(i, 3), "MLG:", length(unique(z))))
    }
  }, img.name = paste0(treefun, "_", measure), imgdir = paste0(measure), 
  htmlfile = paste0(treefun, "_", measure, ".html"), 
  autobrowse = FALSE, title = paste("mlg.filter option", measure), 
  description = desc)
  setwd(cwd)
}
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
set.seed(20150412)
x <- lapply(1:10, getSims, n = 20, snps = 1e4, strucrat = 1, ploidy = 2, err = 0.1, clone = TRUE, n.cores = 4)
y <- lapply(1:10, getSims, n = 20, snps = 1e4, strucrat = 1, ploidy = 2, err = 0.1, clone = FALSE, n.cores = 4)
# x <- getSims(n = 200, snps = 1e4, strucrat = 1, ploidy = 2, err = 0.1, clone = TRUE, n.cores = 4)
# y <- getSims(n = 200, snps = 1e4, strucrat = 1, ploidy = 2, err = 0.1, clone = FALSE, n.cores = 4)
#'
#' For analysis, $\frac{1}{5}$th of the pooled samples will be kept.
x <- do.call("rbind", c(x, y))
x <- x[sample(nInd(x), nInd(x)/5)]
trueclones <- duplicated(substr(indNames(x), start = 1, stop = 10))
fstats <- filter_stats(x, bitwise.dist, plot = TRUE, nclone = sum(!trueclones))

thresh     <- duplicated(mlg.filter(x, distance = bitwise.dist, 
                                    threshold = fstats$farthest$thresholds[sum(trueclones)] + .Machine$double.eps^0.5, 
                                    algo = "f"))
table(thresh, trueclones)
#' 
#' The tabulation is a power analysis of how many true and false positives there
#' are when collapsing at the threshold that gives the same number of known
#' clones/replicates.
#'
#' Below is labelling a tree with known clones.
#+ fig.width = 10, fig.height = 20
the_tree <- upgma(bitwise.dist(x))
clones <- substr(the_tree$tip.label[thresh], start = 1, stop = 10)
clones <- lapply(clones, grep, the_tree$tip.label)
edgelist <- length(which.edge(the_tree, the_tree$tip.label))
edgecols <- rep("black", edgelist)
for (i in clones){
  edgecols[which.edge(the_tree, the_tree$tip.label[i])] <- "red"
}
plot.phylo(the_tree, edge.color = edgecols, adj = 0, label.offset = 0.001)
axisPhylo(1)
title("Random sequences with 10k SNPs and a 0.1 error rate")
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