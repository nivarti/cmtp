# makefile for latex report

all: 	report pdf

full: 	clean all

clean:
	-rm *.ps
	-rm *.pdf
	-rm *.aux
	-rm *.toc
	-rm *.log
	-rm *.dvi
	-rm *.bbl
	-rm *.blg

report:
	latex nse
	latex nse

pdf:
	dvips -o nse.ps nse.dvi
	ps2pdf nse.ps
