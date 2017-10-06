PLT_FILES := $(wildcard *.plt)
PLT_TARGETS := $(PLT_FILES:%.plt=%.tex)
PLT_PDF := $(PLT_FILES:%.plt=%.pdf)
PLT_DVI := $(PLT_FILES:%.plt=%.dvi)
PLT_PS := $(PLT_FILES:%.plt=%.ps)
PLT_EPS := $(PLT_FILES:%.plt=%.eps)

LATEX = latex #-interaction=nonstopmode -file-line-error
DEL = rm

all: $(PLT_PDF) $(PLT_EPS)
.PHONY: clean


FORCE:

%.pdf: %.ps
	ps2pdf $<

%.eps: %.ps
	ps2eps -f $< $@

%.ps: %.dvi %.aux %.log
	dvips $<

%.dvi: %.tex %-inc.eps
	$(LATEX) $<

%.aux:
	@echo ""

%.log:
	@echo ""

%.tex: %.plt
	gnuplot $<

%-inc.eps:
	@echo ""

clean:
	rm -f *.aux
	rm -f *.tex
	rm -f *.dvi
	rm -f *.log
	rm -f *.ps
	rm -f *-inc.eps
	rm -f *eps.eps
