!#/bin/bash
pdflatex test
bibtex test
pdflatex test
evince test.pdf &
