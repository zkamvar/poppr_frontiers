RDIR = .
KNIT = R --slave -e "rmarkdown::render('$<')"

ARTICLE_DIR = $(RDIR)/main_article
RMD_FILES := $(wildcard $(ARTICLE_DIR)/*.Rmd)
TEX_FILES := $(patsubst %.Rmd, %.tex , $(RMD_FILES))
UTEX_FILES := $(patsubst %.tex, %_unicode.tex, $(TEX_FILES))
PDF_FILES := $(patsubst %.tex, %.pdf, $(UTEX_FILES))

maketex: $(TEX_FILES)

%.tex: %.Rmd
	$(KNIT)

%_unicode.tex: %.tex
	python $(ARTICLE_DIR)/convert_pandoc_latex.py '$<' '$@';
	
%.pdf: %.tex
	cd $(ARTICLE_DIR); \
	latexmk -pdf -quiet *_unicode.tex;
	
pdf2eps: 
	sh pdf2eps.sh

unitex: $(UTEX_FILES)
	perl -p -i -e 's/(Figure-\d)-1+?\}/$$1.eps\}/' '$<'
	
difftex: unitex
	cd $(ARTICLE_DIR); \
	latexdiff poppr_frontiers_unicode.tex poppr_frontiers_revision_one_unicode.tex > diff.tex; \
	latexmk -pdf -quiet diff.tex; \
	

final_pdf: $(PDF_FILES)

tex: final_pdf clean

.PHONY: clean

clean: 
	$(RM) $(ARTICLE_DIR)/*.log $(ARTICLE_DIR)/*.out $(ARTICLE_DIR)/*.aux;\
	$(RM) $(ARTICLE_DIR)/*.toc $(ARTICLE_DIR)/*.blg $(ARTICLE_DIR)/*.bbl;\
	$(RM) $(ARTICLE_DIR)/*.synctex.gz $(ARTICLE_DIR)/*.fdb_latexmk $(ARTICLE_DIR)/*.fls;


#===============

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



