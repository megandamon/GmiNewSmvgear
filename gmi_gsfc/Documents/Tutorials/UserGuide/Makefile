
all: 
	latex gmiMain
	bibtex gmiMain
	latex gmiMain
	latex gmiMain
	dvipdf gmiMain.dvi gmiUserGuide.pdf

clean:
	rm -f *.aux *.log *.lof *.dvi *.toc 
	rm -f *.lot *.bbl *.blg *.pdf
