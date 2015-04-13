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