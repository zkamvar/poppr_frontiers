#!/bin/sh

PDFS=main_article/poppr_frontiers_files/figure-latex/*pdf

for i in $PDFS
do
	eps=$(echo $i | rev | cut -b 7- | rev)".eps";
	echo $eps;
	gs -q -dNOCACHE -dNOPAUSE -dBATCH -dSAFER -sDEVICE=epswrite -sOutputFile=$eps $i
done