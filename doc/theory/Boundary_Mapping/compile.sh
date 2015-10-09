#!/bin/sh
pdflatex Boundary_Mapping.tex
bibtex Boundary_Mapping.aux
pdflatex Boundary_Mapping.tex
pdflatex Boundary_Mapping.tex