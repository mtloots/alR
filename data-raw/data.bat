set LaTeXFile=data
rterm --no-save < %LaTeXFile%.r
::pdflatex %LaTeXFile%.tex
::pdflatex %LaTeXFile%.tex
::bibtex %LaTeXFile%
::pdflatex %LaTeXFile%.tex
::pdflatex %LaTeXFile%.tex
::pdflatex %LaTeXFile%.tex

del %LaTeXFile%.aux
del %LaTeXFile%.bbl
del %LaTeXFile%.blg
del %LaTeXFile%.out
del %LaTeXFile%.toc