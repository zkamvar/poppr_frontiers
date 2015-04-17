boot:
	R --slave -e 'install.packages(c("adegenet", "ape", "phangorn", "pegas", "knitr", "animate", "devtools"), repos = "http://cran.at.r-project.org")'; \
	R --slave -e "install.packages('poppr_current.tar.gz', type = 'source', repos = NULL)"; \
	R --slave -e 'devtools::install("FrontiersTemplate/")'
	