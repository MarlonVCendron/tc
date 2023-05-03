filename=documento

all: compile
	
compile:
	latexmk -pdf -pv $(filename).tex -f
	latexmk -c
	@rm $(filename).bbl
	
clean:
	latexmk -c
	@rm $(filename).bbl

watch:
	latexmk -pdf -pvc $(filename).tex