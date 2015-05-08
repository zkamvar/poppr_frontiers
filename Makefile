boot: inst_deps ver_two_deps inst_poppr_two build_template

all: unitex pdf2eps tex

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

build_manuscript:
	cd main_article; \
	R --slave -e 'rmarkdown::render("poppr_frontiers.Rmd")'

uni2tex: build_manuscript
	cd main_article; \
	python convert_pandoc_latex.py poppr_frontiers.tex poppr_frontiers_unicode.tex;

pdf2eps: 
	sh pdf2eps.sh

epstex: uni2tex
	cd main_article; \
	perl -p -i -e 's/(Figure-\d)-1+?\}/$$1.eps\}/' poppr_frontiers_unicode.tex

unitex: epstex

final_pdf:
	cd main_article; \
	latexmk -pdf -quiet poppr_frontiers_unicode.tex;

tex: final_pdf clean

clean: 
	cd main_article; \
	$(RM) *.log *.out *.aux *.toc *.blg *.bbl *.synctex.gz \
	*.fdb_latexmk *.fls
