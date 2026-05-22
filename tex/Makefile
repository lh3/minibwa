.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

all:minibwa.pdf

minibwa.pdf:minibwa.tex minibwa.bib fig-perf.pdf
		pdflatex minibwa; bibtex minibwa; pdflatex minibwa; pdflatex minibwa;

fig-perf.eps:plot-perf.gp perf.tsv
		gnuplot $<

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ minibwa.pdf
