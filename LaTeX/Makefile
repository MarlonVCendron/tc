filename=tc

all: compile

compile:
	latexmk -pdf $(filename).tex
	bibtex $(filename).aux
	latexmk -pdf $(filename).tex
	latexmk -pdf $(filename).tex
	latexmk -c
	@rm $(filename).bbl
	
clean:
	latexmk -c
	@rm $(filename).bbl

watch:
	latexmk -pdf -pvc $(filename).tex