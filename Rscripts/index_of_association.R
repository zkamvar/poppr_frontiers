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
set.seed(999)
x <- glSim(100, 1e3, n.snp.struc=100, ploid=2, alpha = 0.1)
y <- glSim(100, 1.1e3, ploid = 2) # random population
plot(x)
plot(y)
#' 
#' If we use the windowing function for $\bar{r}_d$, here, we can see that all 
#' of the linkage is at the end of the data. First, we are going to set some
#' arbitrary positions for the data.
#+ first_ia, cache = TRUE
set.seed(999)
# I've wanted to do this type of assignment for a long time :)
position(x) <- sort(sample(1.1e4, 1.1e3)) -> position(y)
system.time(x.ia <- win.ia(x, quiet = TRUE)) # window size = 100
system.time(y.ia <- win.ia(y, quiet = TRUE)) # window size = 100
#'
#' Now we can visualize it
plot(x.ia, type = "l", main = "Index of Association over 100nt windows",
     ylab = "Index of Association", xlab = "Window")
lines(y.ia, col = "blue")
legend("topleft", lty = 1, col = c("black", "blue"), legend = c("structured snps", "unstructured snps"))
#'
#' Of course, with this strong LD, it's pretty easy to detect a >100bp chunk.
#' What happens if we randomly shuffle the loci?
set.seed(999)#
shuffled_loci <- sample(nLoc(x))
x.shuff <- x[sample(nInd(x)), shuffled_loci]
position(x.shuff) <- position(x) # Reset the position
#' 
#+ shuffle_ia, cache = TRUE
plot(x.shuff) #
system.time(x.ia.shuff <- win.ia(x.shuff, quiet = TRUE))
#'
#' Now we can plot it!
plot(x.ia.shuff, type = "l", main = "Index of Association over 100nt windows",
     ylab = "Index of Association", xlab = "Window")
abline(h = 0)
#' 
#' We can see exactly where the linked loci are by figuring out what windows they
#' are in:
nwin <- ceiling(max(position(x)/100L))
winmat <- matrix(100L * 1:nwin, nrow = nwin, ncol = 2)
winmat[, 1] <- winmat[, 1] - 100L + 1
linked <- position(x)[shuffled_loci > 1000]
link_wins <- unlist(sapply(linked, 
                           function(x) which(x >= winmat[, 1] & x <= winmat[, 2])))
#' Now we can visualize this all together!
plot(x.ia.shuff, type = "l", main = "Index of Association over 100nt windows",
     ylab = "Index of Association", xlab = "Window")
win_tab <- tabulate(link_wins)
which_wins <- win_tab >= 3
abline(v = which(which_wins), lty = 3, lwd = win_tab[which_wins] - 1)
abline(h = 0)
#'
#'
#' Randomly sampling snps
#' ==========
#' 
#' Another method that can be done is to randomly sample SNPs to see if there
#' is any structure. Here, we are sampling 50 snps 100 times.
#+ snp_samp, cache = TRUE
x.samp <- samp.ia(x, quiet = TRUE, n.snp = 50)
y.samp <- samp.ia(y, quiet = TRUE, n.snp = 50)
#'
#' Now we can produce a boxplot
boxplot(data.frame(structure = x.samp, no_structure = y.samp), 
        ylab = "Index of Association", 
        main = "Index of association over 100 random samples of 50nt")
#'
#' With clones
#' ===========
#' 
#' Now, we can compare populations with clones.
source("my_functions.R")
#'
#+ data_setup_clone, cache = TRUE
set.seed(20150501)
clone_samples <- getSims(n = 100, snps = 1.1e3, strucrat = 1, ploidy = 2, 
                         na.perc = 0.01, err = 0.01, clone = TRUE, n.cores = 4)
set.seed(20150501)
sex_samples <- getSims(n = 100, snps = 1.1e3, strucrat = 1, ploidy = 2, 
                       na.perc = 0.01,  err = 0.01, clone = FALSE, n.cores = 4)
#'
position(clone_samples) <- position(x) -> position(sex_samples)
plot(clone_samples)
plot(sex_samples)
#'
#+ clone_window, cache = TRUE
system.time(clone.ia <- win.ia(clone_samples, quiet = TRUE))
system.time(sex.ia <- win.ia(sex_samples, quiet = TRUE))
#'
plot(clone.ia, type = "l", main = "Index of Association over 100nt windows",
     ylab = "Index of Association", xlab = "Window")
lines(sex.ia, col = "blue")
legend("topleft", lty = 1, col = c("black", "blue"), legend = c("clonal snps", "sexual snps"))
abline(h = mean(clone.ia), lty = 2)
abline(h = mean(sex.ia), lty = 2, col = "blue")
#'
#' Now, to do boxplots
#+ clone_samp, cache = TRUE
set.seed(20150501)
system.time(clone.samp <- samp.ia(clone_samples, quiet = TRUE, n.snp = 50))
set.seed(20150501)
system.time(sex.samp <- samp.ia(sex_samples, quiet = TRUE, n.snp = 50))
#' Now we can produce a boxplot
boxplot(data.frame(clonal = clone.samp, sexual = sex.samp), 
        ylab = "Index of Association", 
        main = "Index of association over 100 random samples of 50nt")
#' ## Session Info
options(width = 100)
devtools::session_info()