boot: inst_deps inst_poppr build_template

inst_deps:
	R --slave -e 'install.packages(c("adegenet", "ape", "phangorn", "pegas", "knitr", "rmarkdown", "animation", "devtools"), repos = "http://cran.at.r-project.org")'; \

inst_poppr:
	R --slave -e "install.packages('poppr_current.tar.gz', type = 'source', repos = NULL)"; \

build_template:
	R --slave -e 'devtools::install("FrontiersTemplate/")'
