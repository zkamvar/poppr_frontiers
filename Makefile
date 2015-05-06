boot: inst_deps ver_two_deps inst_poppr_two build_template

inst_deps:
	R --slave -e 'install.packages(c("ape", "knitr", "rmarkdown", "animation", "devtools"), repos = "http://cran.at.r-project.org")'; \

inst_poppr_one:
	R --slave -e "devtools::install_github('grunwaldlab/poppr')"; \

inst_poppr_two:
	R --slave -e "devtools::install_github('grunwaldlab/poppr@2.0-rc')";

build_template:
	R --slave -e 'devtools::install("FrontiersTemplate/")'

ver_two_deps:
	R --slave -e 'devtools::install_github(c("thibautjombart/adegenet", "emmanuelparadis/pegas/pegas", "KlausVigo/phangorn"))';

unitex:
	cd main_article; \
	python convert_pandoc_latex.py poppr_frontiers.tex poppr_frontiers_unicode.tex;

tex:
	cd main_article; \
	latexmk -pdf -quiet poppr_frontiers_unicode.tex;

clean: 
	cd main_article; \
	$(RM) *.log *.out *.aux *.toc *.blg *.bbl *.synctex.gz \
	*.fdb_latexmk *.fls
