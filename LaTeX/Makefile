filename=tc

all: compile

compile:
	latexmk -pdf -pv $(filename).tex
	bibtex $(filename).aux
	latexmk -pdf -pv $(filename).tex
	latexmk -pdf -pv $(filename).tex
	latexmk -c
	@rm $(filename).bbl
	
clean:
	latexmk -c
	@rm $(filename).bbl

watch:
	latexmk -pdf -pvc $(filename).tex