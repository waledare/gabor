!#/bin/bash
pdflatex firstpresjp
bibtex firstpresjp
pdflatex firstpresjp
evince firstpresjp.pdf &
