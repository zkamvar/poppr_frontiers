# Frontiers Article for the updated poppr


# INSTRUCTIONS FOR SETUP (**updated**)

1. Open `poppr_frontiers.Rproj` in Rstudio.
2. Run `make boot`. This will do the following:
    i. Install dependent packages.
    ii. Install the version of poppr necessary for this.
    iii. install the FrontiersTemplate.
4. Knit the file `main_article/poppr_frontiers.Rmd`

## To prepare file document

1. run `make unitex`. This will convert the unicode tex document to pure tex and then convert pdf figures to eps and finally replace the figure names with eps names in the tex document.

2. move tables and figures to the bottom of the manuscript.

3. run `make difftex`. This will create a diff file between the original and revision.

4. run `make final_pdf` This will make the final latex document called poppr_frontiers_unicode.pdf

# Files and Directories

 - Frontiers\_LaTeX\_Templates - original template from Frontiers
 - FrontiersTemplate - Rmarkdown template based off of http://github.com/rstudio/rticles
 - images - where visualizations are stored for the time being
 - main_article - where the paper is written
 - Rscripts - various R scripts

********

This repository will contain the writeup for the frontiers article for the improved version of poppr that includes Jonah Brooks as an author.

## What this article needs to highlight

- The rise in popularity of reproducible methods for population genetic research.
    - moving away from standalone programs and into R.
- Population genetics packages in R
- The lack of tools for analysis of clonal populations and population hierarchies
    - poppr formally implements hierarchies **moved to adegenet, but can still talk about**
- Improvements to poppr since its inception

## Poppr's improvements

Note that when we talk about the improvements to poppr, we are talking about the improvements after version 1.0.

- Formal methods to deal with hierarchical levels.
- Multilocus genotype flexibility (different word, maybe)
    - Collapsing by genetic distance (addressing issues with sequencing error)
    - Expanding by Psex value (addressing issues with too few markers) **Future**
    - user-defined version **Future**
- Psex
- Bootstrapping of any genetic distance measure for individuals and populations. 
- Index of association for genomic data (windowing, random sampling of loci, whole-genome)
- Parallelization
- Genotype accumulation curve
- implementation of loops in MSN for analysis with graph-walking algorithms (in igraph)