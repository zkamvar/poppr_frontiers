getNames <- function(n){
  vapply(1:n, function(x) paste(sample(letters, 10), collapse = ""), character(1))
}

# This will mutate a single allele in a random chromosome when given a SNPbin
snp_mutator <- function(chrom, mutation = sample(2^(1:7), 1)){
  posi        <- sample(length(chrom) - 1, 1) # sample chunk of 8 loci
  orig        <- as.integer(chrom[posi])      # convert to integer
  chrom[posi] <- as.raw(bitwXor(orig, mutation)) # Exclusive OR will flip the bit. 
  return(chrom)
}

# This will mutate nmutations in a SNPbin
sample_mutator <- function(snpbin, mu, nLoc, rawchars = 2^(0:7)){
  nmutations <- rpois(1, lambda = round(nLoc*mu))
  for (i in seq(nmutations)){
    chrom_index               <- sample(length(snpbin@snp), 1)
    chrom                     <- snpbin@snp[[chrom_index]]
    snpbin@snp[[chrom_index]] <- snp_mutator(chrom, sample(rawchars, 1))
  }
  return(snpbin)
}

# This will mutate 
pop_mutator <- function(glt, mu = 0.05, samples = TRUE){
  rawchrs <- 2^(0:7)
  glt@gen[samples] <- lapply(glt@gen[samples], sample_mutator, mu, nLoc(glt), rawchrs)
  return(glt)
}

NA_zeromancer <- function(chrom, NA.posi, rawchars = 2^(0:7)){
  nas <- ceiling(NA.posi/8) # Getting the location in the RAW vector
  zero_bits <- NA.posi %% 8 # Getting the location of the locus in a RAW element.
  zero_bits[zero_bits == 0] <- 8
  
  for (i in seq(length(nas))){
    eight_bit <- as.integer(chrom[nas[i]])
    the_locus <- rawchars[zero_bits[i]]
    # If the locus does not change in the OR, then there is a 1 and it needs to
    # be changed to a zero with XOR.
    if (eight_bit == bitwOr(eight_bit, the_locus)){
      chrom[nas[i]] <- as.raw(bitwXor(eight_bit, the_locus))
    }
  }
  return(chrom)
}

NA_generator <- function(snpbin, nloc, na.perc = 0.01, rawchars = 2^(0:7)){
  nas <- rpois(1, lambda = round(nloc*na.perc))
  NA.posi <- sort(sample(nloc, nas))
  for (i in seq(length(snpbin@snp))){
    snpbin@snp[[i]] <- NA_zeromancer(snpbin@snp[[i]], NA.posi, rawchars)
  }
  snpbin@NA.posi <- NA.posi
  return(snpbin)
}

pop_NA <- function(glt, na.perc = 0.01, parallel = require('parallel'), n.cores = 2L){
  rawchars <- 2^(0:7)
  if (parallel){
    glt@gen <- mclapply(glt@gen, NA_generator, nLoc(glt), na.perc, rawchars,
                        mc.cores = getOption("mc.cores", n.cores))
  } else {
    glt@gen <- lapply(glt@gen, NA_generator, nLoc(glt), na.perc, rawchars)    
  }

  return(glt)
}

crossover <- function(snpbin){
  chr1 <- snpbin@snp[[1]]
  chr2 <- snpbin@snp[[2]]
  chrlen <- length(chr1)
  cutpoint <- sample(2:(chrlen - 1), 1)
  first <- 1:(cutpoint - 1)
  last  <- cutpoint:chrlen
  snpbin@snp[[1]] <- c(chr1[first], chr2[last])
  snpbin@snp[[2]] <- c(chr2[first], chr1[last])
  return(snpbin)
}

mate <- function(snpbin, ind1, ind2){
  snpbin@gen[[ind1]] <- crossover(snpbin@gen[[ind1]])
  snpbin@gen[[ind2]] <- crossover(snpbin@gen[[ind2]])
  snpout <- snpbin@gen[[ind1]]
  snpout@snp[[1]] <- snpbin@gen[[ind1]]@snp[[sample(1:2, 1)]]
  snpout@snp[[2]] <- snpbin@gen[[ind2]]@snp[[sample(1:2, 1)]]
  return(snpout)
}

random_mate <- function(glt, err){
  mat_pair <- matrix(integer(1), nrow = nInd(glt), ncol = 2)
  mat_pair[, 1] <- sample(nInd(glt), replace = TRUE)
  mat_pair[, 2] <- sample(nInd(glt), replace = TRUE)
  res <- apply(mat_pair, 1, function(x) mate(glt, x[1], x[2]))
  glt@gen <- res
  if (err > 0) glt <- pop_mutator(glt, err)
  return(glt)
}

random_mate_gen <- function(glt, err = 5e-3, gen = 1){
  for (i in seq(gen)){
    glt <- random_mate(glt, err)
  }
  return(glt)
}


getSims <- function(z = 1, n = 10, snps = 1e6, strucrat = c(0.25, 0.75), 
                    clone = TRUE, err = 0.1, na.perc = 0.1, n.cores = 4, 
                    mate_gen = NULL, mate_err = 5e-3, ...){
  lam <- sample(strucrat, 1)
  if (lam == 1){
    perc <- 1
  } else {
    perc <- ceiling(sample(rpois(1000, lambda = snps*lam), 1)/snps)    
  }

  n <- sample(rpois(1000, lambda = n), 1)
  the_names <- getNames(n)
  res <- glSim(n, n.snp.nonstruc = snps*perc, n.snp.struc = snps*(1-perc), 
               n.cores = n.cores, ...)
  if (!is.null(mate_gen)){
    res <- random_mate_gen(res, mate_err, mate_gen)
  }
  if (clone){
    clones <- sample(n, replace = TRUE)
    the_names <- paste(the_names[clones], 1:n, sep = ".")
    res <- res[clones]
    clones <- duplicated(clones)
    # res <- pop_mutator(res, err, clones)
  }
  if (err > 0){
    res <- pop_mutator(res, err)    
  }
  if (na.perc > 0){
    res <- pop_NA(res, na.perc = na.perc, n.cores = n.cores)    
  }
  indNames(res) <- the_names
  return(res)
}

binary_char_from_hex <- function(y){
  vapply(y, function(x) paste(as.integer(rawToBits(x)), collapse=""), character(1))
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


threshold_predictor <- function(thresholds, fraction = 0.5){
  frac <- 1:round(length(thresholds)*fraction)
  diffs <- diff(thresholds[frac])
  diffmax <- which.max(diffs)
  mean(thresholds[diffmax:(diffmax + 1)])
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