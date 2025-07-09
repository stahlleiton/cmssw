pdf : report.pdf
dvi : report.dvi

INTER = -interaction batchmode

%.pdf : %.tex $(wildcard ../eps/*.eps) $(wildcard ../eps/*/*.eps)
	pdflatex $(INTER) $<
	pdflatex $(INTER) $<

%.dvi : %.tex $(wildcard ../eps/*.eps) $(wildcard ../eps/*/*.eps)
	latex $(INTER) $<
	latex $(INTER) $<
	killall --signal SIGUSR1 xdvi.bin &

clean :
	rm -f report.dvi report.pdf report.aux report.log report.out report.toc
