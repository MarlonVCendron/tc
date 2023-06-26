filename=tc

all: compile

compile:
	latexmk -pdf $(filename).tex
	bibtex $(filename).aux
	latexmk -pdf $(filename).tex
	latexmk -pdf $(filename).tex
	latexmk -c
	@rm -f $(filename).bbl
	@rm -f $(filename).brf
	
clean:
	latexmk -c
	@rm -f $(filename).bbl
	@rm -f $(filename).brf

watch:
	latexmk -pdf -pvc $(filename).tex