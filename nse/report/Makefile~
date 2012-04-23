
# The Latex fruit poster template
# http://www.dcs.qmul.ac.uk/~andrew/poster/
# andrew graves 2004 andrew@dcs.qmul.ac.uk

all: 	poster pdf

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
	-rm html/*
	-rmdir html

poster:
	latex poster
	latex poster

pdf:
	dvips -o poster.ps poster.dvi
	ps2pdf -dAutoFilterColorImages=false -sColorImageFilter=FlateEncode poster.ps

jpg:
	convert -geometry 25% poster.pdf poster.jpg

