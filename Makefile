FIGURES=\
	figures/tree-sequence-illustration.pdf\
	figures/tree-transition.pdf\
	figures/otex.pdf\
	figures/num_events.pdf	

paper.pdf: paper.tex paper.bib ${FIGURES}
	pdflatex paper.tex
	bibtex paper
	pdflatex paper.tex

paper.ps: paper.dvi 
	dvips paper

paper.dvi: paper.tex paper.bib 
	latex paper.tex
	bibtex paper
	latex paper.tex
	latex paper.tex

figures/tree-sequence-illustration.pdf: src/illustrations.py
	python src/illustrations.py

figures/otex.pdf: src/otex.asy 
	asy -o figures/otex.pdf -fpdf src/otex.asy
	asy -o figures/otex.eps -feps src/otex.asy

figures/tree-transition.pdf: src/tree-transition.asy 
	asy -o figures/tree-transition.pdf -fpdf src/tree-transition.asy
	asy -o figures/tree-transition.eps -feps src/tree-transition.asy

figures/num_events.pdf: src/plots.py
	python src/plots.py plot all

clean:
	rm -f *.log *.dvi *.aux
	rm -f *.blg *.bbl
	rm -f figures/*
	rm -f *.eps 
	rm -f tags

tags: paper.tex paper.bib src/*.py
	ctags src/*.py # paper.tex paper.bib 


mrproper:
	rm -f *.ps *.pdf
