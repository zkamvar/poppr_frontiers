boot: inst_deps inst_ver_two build_template

inst_deps:
	R --slave -e 'install.packages(c("ape", "knitr", "rmarkdown", "animation", "devtools"), repos = "http://cran.at.r-project.org")'; \

inst_poppr:
	R --slave -e "devtools::install_github('grunwaldlab/poppr')"; \

build_template:
	R --slave -e 'devtools::install("FrontiersTemplate/")'

inst_ver_two:
	R --slave -e 'devtools::install_github(c("thibautjombart/adegenet", "emmanuelparadis/pegas/pegas", "KlausVigo/phangorn", "grunwaldlab/poppr@adegenet-fix"))'; \
