!#/bin/bash
pdflatex first
bibtex first
pdflatex first
evince first.pdf &
