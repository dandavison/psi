
BASE = emvbpl
AUX = $(BASE).aux
BBL = $(BASE).bbl
TEX = $(BASE).tex
JUNK = $(BASE).toc $(BASE).out $(BASE).log $(BASE).blg $(BASE).dvi

LATEX = pdflatex
BIBTEX = bibtex


$(AUX):
	$(LATEX) $(BASE) > /dev/null

$(BIB): $(AUX)
	bibtex $(BASE)

pdf: 	$(BIB)
	make clean
	$(LATEX) $(BASE) > /dev/null
	bibtex $(BASE)
	$(LATEX) $(BASE) > /dev/null
	$(LATEX) $(BASE) > /dev/null

clean:
	rm -f $(AUX) $(BBL) $(JUNK)
