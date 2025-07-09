EPSFILES = $(wildcard *.eps) $(wildcard */*.eps)
PDFFILES = $(patsubst %.eps,%.pdf,$(EPSFILES))

all : $(PDFFILES)

%.pdf : %.eps
	ps2pdf -dEPSCrop $< $@

clean :
	rm -f *.eps */*.eps
	rm -f *.pdf */*.pdf
