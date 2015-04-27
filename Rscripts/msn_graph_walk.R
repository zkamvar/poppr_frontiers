#' Demonstrating Graph Walking algorithms in R
#' =============
#' 
#' In version 2.0, poppr introduces the option to include reticulation in 
#' minimum spanning networks. This is essential for clonal organisms as mutation
#' has a much higher influence on genetic diversity as opposed to meiotic 
#' processes. 
#' 
#' In this examle, I will utilize *P. ramorum* data from Kamvar et al. (2015). The source
#' should be downloaded with the following command:
#+ eval = FALSE
devtools::install_github("zkamvar/Sudden_Oak_Death_in_Oregon_Forests/PramCurry")
#' 
#' After you've downloaded the package, you can load the library and the data.
library("PramCurry")
library("poppr")
library("igraph")
data(for2nur)
data(comparePal)
#' 
#' Below are some functions we need to convert the data to the proper format as
#' well as some helper functions from the PramCurry package.
# Conversion from adegenet 1.4 to 2.0
data_replacer <- function(x){
  # manipulate your data here
  if (x@type == "codom"){
    xtab <- x@tab
    newtab <- as.integer(x@tab * x@ploidy)
    x@tab <- matrix(newtab, nrow = nrow(xtab), ncol = ncol(xtab),
                    dimnames = dimnames(xtab))
    x@ploidy <- rep(x@ploidy, nrow(x@tab))
    names(x@loc.names) <- NULL
    rownames(x@tab) <- x@ind.names
    names(x@all.names) <- x@loc.names
    colnames(x@tab) <- unlist(lapply(x@loc.names, 
                                     function(i) paste(i, x@all.names[[i]], 
                                                       sep = ".")), 
                              use.names = FALSE)
    levels(x@loc.fac) <- x@loc.names
    names(x@loc.nall) <- x@loc.names
  }
  
  if (!is.null(x@pop)){
    levels(x@pop) <- x@pop.names
  }
  if ("genclone" %in% class(x)){
    x@strata <- x@hierarchy
    x@hierarchy <- NULL
  } else {
    x@hierarchy <- NULL
    x@strata    <- NULL
  }
  return(x)
}

# These functions are from Kamvar et al. 2015
make_node_list <- function(clusts, pal, cutoff = 3){
  PAL <- match.fun(pal)
  clustPal  <- table(clusts$membership)
  theClusts <- clustPal > 3
  clustOrder <- order(clustPal, decreasing = TRUE)
  clustPal[clustOrder][theClusts[clustOrder]]   <- PAL(sum(theClusts))
  clustPal[clustOrder][!theClusts[clustOrder]]  <- gray.colors(sum(!theClusts))
  nodeList <- lapply(1:length(clustPal), function(x) which(clusts$membership == x))
  names(nodeList) <- clustPal
  return(nodeList)
}

clust_cutoff <- function(clusts, cutoff = 3){
  table(clusts$membership) > 3
}
#' 
#' ### Data setup
#' Conversion
for2nur <- data_replacer(for2nur)
setPop(for2nur) <- ~SOURCE/STATE
for2nur
#' 
#' Repeat lengths are necessary:
newReps <- other(for2nur)$REPLEN
(newReps[3] <- 4) # Tetranucleotide repeat
(newReps <- fix_replen(for2nur, newReps))
#'
#' ### Calculating the minimum spanning networks
#' 
#' We will create two msns: one with no reticulations (`nf.msn.noties`) and
#' one with reticulations (`nf.msn.ties`). these will be analyzed with the
#' infomap community detection algorithm and then the memberships will be compared
#' to the observed populations. 
# No reticulations in the network.
nf.msn.noties <- bruvo.msn(for2nur, replen = newReps, showplot = FALSE)
# Adding reticulations
nf.msn.ties <- bruvo.msn(for2nur, replen = newReps, include.ties = TRUE,
                         showplot = FALSE)
#' 
#' Now, we are going to use the igraph function `infomap.community`. Basically,
#' it will send random walkers through the graph and see where they spend the
#' most time. It considers edge weights and vertex weights. For vertex weights,
#' we are setting them to be the number of samples in each vertex. 
clusts.noties <- infomap.community(nf.msn.noties$graph, nb.trials = 1e3, 
                                   v.weights = V(nf.msn.noties$graph)$size)
clusts.ties <- infomap.community(nf.msn.ties$graph, nb.trials = 1e3, 
                                 v.weights = V(nf.msn.ties$graph)$size)
#'
#' The membership of each sample is matched here. 
members.ties <- clusts.ties$membership[match(mlg.vector(for2nur), mlgFromString(V(nf.msn.ties$graph)$label))]
members.noties <- clusts.noties$membership[match(mlg.vector(for2nur), mlgFromString(V(nf.msn.noties$graph)$label))]
#'
#' Now, we are comparing the memberships vs. the populations
table.value(table(members.ties, pop(for2nur)), col.labels = popNames(for2nur))
table.value(table(members.noties, pop(for2nur)), col.labels = popNames(for2nur))
#'
#' One other thing to look at is the sizes of the communities. With a clonal
#' pathogen, we expect to see larger communities. 
sizes(clusts.ties)
sizes(clusts.noties)
#'
#' What we are seeing is that the reticulations allow for more genotypes to be
#' clustered together because there are more meaningful connections. 
#' 
#' Below, we will see what the graphs look like with and without ties.
#' ### with ties
#+ ties, fig.width = 10, fig.height = 10
nf.msn <- nf.msn.ties
clusts <- clusts.ties
mypal <- ifelse(length(unique(clusts$membership)) < 14, 
                RColorBrewer::brewer.pal(x, "Paired"), 
                colorRampPalette("black"))
nodeList <- make_node_list(clusts, mypal, cutoff = 3)
theClusts <- clust_cutoff(clusts, 3)

goodSeed <- 6
thisPal <- function(x) comparePal[nf.msn$populations]

set.seed(goodSeed)
MASTER <- get_layout(nf.msn$graph, LAYOUT = layout.fruchterman.reingold)
plot_poppr_msn(for2nur, nf.msn, gad = 10, palette = thisPal, mlg = TRUE, 
               layfun = MASTER, nodebase = 1.75, vertex.label.font = 2,
               quantiles = FALSE, #inds = "none",
               vertex.label.color = "firebrick",
               mark.groups = nodeList[theClusts], 
               mark.border = names(nodeList)[theClusts],
               mark.col = transp(names(nodeList)[theClusts], 0.05),
               mark.expand = 2,
               mark.shape = 0)
plot(clusts, nf.msn$graph, vertex.size = log(V(nf.msn$graph)$size, 1.75) + 3, 
     vertex.label = NA, layout = MASTER)

#' ### without ties
#+ noties, fig.width = 10, fig.height = 10
nf.msn <- nf.msn.noties
clusts <- clusts.noties
mypal <- ifelse(length(unique(clusts$membership)) < 14, 
                RColorBrewer::brewer.pal(x, "Paired"), 
                colorRampPalette("black"))
nodeList <- make_node_list(clusts, mypal, cutoff = 3)
theClusts <- clust_cutoff(clusts, 3)

goodSeed <- 6
thisPal <- function(x) comparePal[nf.msn$populations]
set.seed(goodSeed)
MASTER <- get_layout(nf.msn$graph, LAYOUT = layout.auto)
plot_poppr_msn(for2nur, nf.msn, gad = 10, palette = thisPal, mlg = TRUE, 
               layfun = MASTER, nodebase = 1.75, vertex.label.font = 2,
               quantiles = FALSE, #inds = "none",
               vertex.label.color = "firebrick",
               mark.groups = nodeList[theClusts], 
               mark.border = names(nodeList)[theClusts],
               mark.col = transp(names(nodeList)[theClusts], 0.05),
               mark.expand = 2,
               mark.shape = 0)
plot(clusts, nf.msn$graph, vertex.size = log(V(nf.msn$graph)$size, 1.75) + 3, 
     vertex.label = NA, layout = MASTER)