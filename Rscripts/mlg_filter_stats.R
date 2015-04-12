getNames <- function(n){
  vapply(1:n, function(x) paste(sample(letters, 10), collapse = ""), character(1))
}

getSims <- function(z = 1, n = 10, snps = 1e6, strucrat = c(0.25, 0.75), 
                    clone = TRUE, err = 0.1, ...){
  lam <- sample(strucrat, 1)
  if (lam == 1){
    perc <- 1
  } else {
    perc <- ceiling(sample(rpois(1000, lambda = snps*lam), 1)/snps)    
  }

  n <- sample(rpois(1000, lambda = n), 1)
  the_names <- getNames(n)
  res <- glSim(n, n.snp.nonstruc = snps*perc, n.snp.struc = snps*(1-perc), ...)
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

  res <- new("genlight", mat, ploidy = 2, ind.names = the_names)
  return(res)
}


filter_stats <- function(x, distance = bitwise.dist, threshold = 1, 
                         stats = "All", missing = "ignore", plot = FALSE, ...){
  f <- mlg.filter(x, threshold, missing, algorithm = "f", distance = distance, stats = stats, ...)
  a <- mlg.filter(x, threshold, missing, algorithm = "a", distance = distance, stats = stats, ...)
  n <- mlg.filter(x, threshold, missing, algorithm = "n", distance = distance, stats = stats, ...)
  fanlist <- list(farthest = f, average = a, nearest = n)
  if (stats == "All"){
    for (i in names(fanlist)){
      names(fanlist[[i]]) <- c("MLG", "thresholds", "mat", "size")
    }
    if (plot){
      xdist <- distance(x, ...)
      upper <- round(max(xdist), digits = 1)
      plot(x = c(upper, 0), y = c(ifelse(is.genind(x), mlg(x, quiet = TRUE), nInd(x)), 1), type = "n",
           ylab = "Number of Multilocus Lineages",
           xlab = "Genetic Distance Cutoff")
      points(x = rev(a[[2]]), y = 1:length(a[[2]]))
      points(x = rev(f[[2]]), y = 1:length(f[[2]]), col = "blue")
      points(x = rev(n[[2]]), y = 1:length(n[[2]]), col = "red")    
    }
  }
  return(fanlist)
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

RAD_data_dir <- "~/Documents/Grunwald/short-scripts/Genotype_error/"

library('poppr')
library('ape')
library('phangorn')
library('animation')
infdat <- RCurl::getURL("https://raw.githubusercontent.com/grunwaldlab/phytophthora_id/master/shiny-server/www/genoid_infestans/reduced_database.txt.csv")
infdat <- read.table(text = infdat, head = TRUE)
pinf <- df2genind(infdat[-c(1,2)], sep = "/", ploidy = 3, ind.names = infdat[[1]], pop = infdat[[2]])
ssr <- c(3,3,2,3,3,2,2,3,3,3,3,3)
x <- as.genclone(pinf)
fstats <- filter_stats(x, bruvo.dist, plot = TRUE, replen = ssr, loss = FALSE)
title(main = expression(paste(italic("P. infestans"), " reference isolates (12 SSR loci)")))
legend("topright", legend = c("Nearest Neighbor", "UPGMA", "Farthest Neighbor"), 
       col = c("red", "black", "blue"), pch = 1, title = "Clustering Method")
z <- filter_stats(x, bruvo.dist, replen = ssr, loss = FALSE, threshold = 0.75/12, stats = "MLGS")
print(table(pop(x), z$farthest), zero.print = ".")
print(table(pop(x), z$nearest), zero.print = ".")
tab <- mlg.table(x, bar = FALSE)
colnames(tab) <- 1:ncol(tab)
print.table(tab, zero.print = ".")



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


for (i in names(fstats)){
  HTML_collapse(measure = i, x = x, treefun = "upgma", 
                distfun = "bruvo.dist", fstats = fstats, destdir = "images", 
                replen = ssr, loss = FALSE)  
  HTML_collapse(measure = i, x = x, treefun = "nj", 
                distfun = "bruvo.dist", fstats = fstats, destdir = "images", 
                replen = ssr, loss = FALSE)  
}

# 
# 
# desc <- "An exploration of collapsing MLGs using UPGMA"
# saveHTML({
#   for (i in unique(fstats$average$thresholds)){
#     z <- filter_stats(x, bruvo.dist, replen = ssr, loss = FALSE, threshold = i, stats = "MLGS")
#     color_mlg_tree(x, pitree, z$average)
#     add.scale.bar()
#     title(paste("Threshold:", round(i, 3), "MLG:", length(unique(z$average))))
#   }
# }, img.name = "UPGMA", imgdir = "upgma", htmlfile = "UPGMA.html", 
# autobrowse = FALSE, title = "mlg.filter option average", 
# description = desc)
# desc <- "An exploration of collapsing MLGs using Nearest Neighbor"
# saveHTML({
#   for (i in unique(fstats$nearest$thresholds)){
#     z <- filter_stats(x, bruvo.dist, replen = ssr, loss = FALSE, threshold = i, stats = "MLGS")
#     color_mlg_tree(x, pitree, z$nearest)
#     add.scale.bar()
#     title(paste("Threshold:", round(i, 3), "MLG:", length(unique(z$nearest))))
#   }
# }, img.name = "NN", imgdir = "nn", htmlfile = "NN.html", 
# autobrowse = FALSE, title = "mlg.filter option nearest", 
# description = desc)
# desc <- "An exploration of collapsing MLGs using Farthest Neighbor"
# saveHTML({
#   for (i in unique(fstats$farthest$thresholds)){
#     z <- filter_stats(x, bruvo.dist, replen = ssr, loss = FALSE, threshold = i, stats = "MLGS")
#     color_mlg_tree(x, pitree, z$farthest)
#     add.scale.bar()
#     title(paste("Threshold:", round(i, 3), "MLG:", length(unique(z$farthest))))
#   }
# }, img.name = "FN", imgdir = "fn", htmlfile = "FN.html", 
# autobrowse = FALSE, title = "mlg.filter option farthest", 
# description = desc)


x <- lapply(1:10, getSims, snps = 1e4, strucrat = 1, ploidy = 2, err = 0.1, clone = TRUE)
y <- lapply(1:10, getSims, snps = 1e4, strucrat = 1, ploidy = 2, err = 0.1, clone = FALSE)
x <- do.call("rbind", c(x, y))
x <- x[sample(nInd(x), nInd(x)/2)]
fstats <- filter_stats(x, bitwise.dist, plot = TRUE)


trueclones <- duplicated(substr(indNames(x), start = 1, stop = 10))
thresh     <- duplicated(mlg.filter(x, distance = bitwise.dist, 
                                    threshold = 0.2, algo = "a"))
table(thresh, trueclones)

plinklist <- list(m3 = character(0), m4 = character(0), m10 = character(0), def = character(0))
plinklist[["m4"]] <- "2R/PopSamples/data.out/PopSamples_m4/Popsouts_Rselec/out.replicates/plink.raw"

plinklist[["def"]] <- "2R/PopSamples/data.out/PopSamples_def/Popsouts_Rselec/out.replicates/plink.raw"

plinklist[["m3"]] <- "2R/PopSamples/data.out/PopSamples_m3/Popsouts_Rselec/out.replicates/plink.raw"

plinklist[["m10"]] <- "2R/PopSamples/data.out/PopSamples_m10/Popsouts_Rselec/out.replicates/plink.raw"

par(mfrow = c(2, 2))

contlist <- plinklist
threshlist <- plinklist

for (i in names(plinklist)){
  barb <- read.PLINK(paste(RAD_data_dir, plinklist[[i]], sep = "/"))
  show(barb)
  fstats <- filter_stats(barb, bitwise.dist, plot = TRUE)
  title(paste(i, "SNPS:", nLoc(barb)))
  minthresh  <- fstats$average$thresholds[10]
  abline(v = minthresh, lty = 2)
  thresh     <- mlg.filter(barb, distance = bitwise.dist, algorithm = "a", 
                           threshold = minthresh)
  trueclones <- vapply(strsplit(indNames(barb), "_"), "[[", character(1), 1)
  trueclones <- duplicated(trueclones)
  thresh     <- duplicated(thresh)
  contlist[[i]] <- table(thresh, trueclones)
  threshlist[[i]] <- fstats$average$thresholds
}
print(contlist)
for (i in threshlist){
  plot(diff(i), log = "y")
}
par(mfrow = c(1, 1))

defupgma <- phangorn::upgma(bitwise.dist(barb))
z <- filter_stats(barb, bitwise.dist, threshold = minthresh, stats = "MLGS")
color_mlg_tree(barb, defupgma, z$farthest)
axisPhylo(1)
