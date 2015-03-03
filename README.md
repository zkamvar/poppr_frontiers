# Frontiers Article for the updated poppr

This repository will contain the writeup for the frontiers article for the improved version of poppr that includes Jonah Brooks as an author.

## What this article needs to highlight

- The rise in popularity of reproducible methods for population genetic research.
    - moving away from standalone programs and into R.
- Population genetics packages in R
- The lack of tools for analysis of clonal populations and population hierarchies
- Improvements to poppr since its inception

## Poppr's improvements

Note that when we talk about the improvements to poppr, we are talking about the improvments after version 1.0.

- Formal methods to deal with hierarchical levels.
- Multilocus genotype flexibility
    - Collapsing by genetic distance (addressing issues with sequencing error)
    - Expanding by psex value (addressing issues with too few markers)
    - user-defined version
- Bootstrapping of any genetic distance measure for individuals and populations. 
- Index of association for genomic data (windowing, random sampling of loci, whole-genome)
- Parallelization

