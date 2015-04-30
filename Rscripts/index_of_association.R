#' Sliding window Index of Association
#' ============
#' 
#' The index of assocation, developed by Brown, refined by Smith, and then 
#' standardized by Agapow and Burt, is a useful tool for detecting multilocus 
#' linkage disequilibrium, which can serve as an indication of clonal
#' reproduction. Developed for a small number of markers, handling $\bar{r}_d$
#' for SNP data, which could have thousands of loci was not clear. Here, we will
#' show calculation of the index of association utilizing sliding windows,
#' detecting linkage between two or more loci per window.
#' 
library("poppr")
#'
#' The data. This comes from the adegenet package, a simulation of 100
#' individuals with 1,000 unstructured loci and 100 structured loci in strong
#' LD.
#+ data_setup, cache = TRUE
set.seed(2015)
x <- glSim(100, 1e3, n.snp.struc=100, ploid=2, alpha=0.4, LD=TRUE, block.minsize=100)
plot(x)
#' 
#' If we use the windowing function for $\bar{r}_d$, here, we can see that all 
#' of the linkage is at the end of the data. First, we are going to set some
#' arbitrary positions for the data.
#+ first_ia, cache = TRUE
set.seed(2015)
position(x) <- sort(sample(1.1e4, 1.1e3))
system.time(x.ia <- win.ia(x, quiet = TRUE)) # window size = 100
#'
#' Now we can visualize it
plot(x.ia, type = "l", main = "Index of Association over 100nt windows",
     ylab = "Index of Association", xlab = "Window")
#'
#' Of course, with this strong LD, it's pretty easy to detect a >100bp chunk.
#' What happens if we randomly shuffle the loci?
set.seed(2015)
shuffled_loci <- sample(nLoc(x))
x.shuff <- x[, shuffled_loci]
position(x.shuff) <- position(x) # Reset the position
#' 
#+ shuffle_ia, cache = TRUE
plot(x.shuff) # 
system.time(x.ia.shuff <- win.ia(x.shuff, quiet = TRUE))
#'
#' Now we can plot it!
plot(x.ia.shuff, type = "l", main = "Index of Association over 100nt windows",
     ylab = "Index of Association", xlab = "Window")
#' 
#' We can see exactly where the linked loci are by figuring out what windows they
#' are in:
nwin <- ceiling(max(position(x)/100L))
winmat <- matrix(100L * 1:nwin, nrow = nwin, ncol = 2)
winmat[, 1] <- winmat[, 1] - 100L
linked <- position(x)[shuffled_loci > 1000]
link_wins <- unlist(sapply(linked, 
                           function(x) which(x >= winmat[, 1] & x <= winmat[, 2])))
two_wins <- table(link_wins)
two_wins <- as.numeric(names(two_wins)[two_wins >= 2])
#' Now we can visualize this all together!
plot(x.ia.shuff, type = "l", main = "Index of Association over 100nt windows",
     ylab = "Index of Association", xlab = "Window")
abline(v = two_wins, lty = 3)
