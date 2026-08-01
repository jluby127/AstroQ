# Build the HIRES 2026B paper.
#
#   make          build main.pdf
#   make params   regenerate params.tex from data/
#   make check    fail if any \param{} key is unresolved in the PDF
#   make clean    remove latexmk build products

PYTHON ?= python
LATEXMK ?= latexmk

.PHONY: all params check clean

all: main.pdf

# params.tex is generated, never hand-edited, so it is rebuilt whenever any
# recorded artifact changes.
params: params.tex

params.tex: tools/make_params.py $(wildcard data/*/*)
	$(PYTHON) tools/make_params.py

main.pdf: main.tex params.tex variablesTable.tex bib.bib
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error main.tex

# An unresolved \param{} key renders as a literal red "XX" in the PDF. Catch it
# at the source level instead, which needs no pdftotext and names the key.
check:
	$(PYTHON) tools/make_params.py --check

clean:
	$(LATEXMK) -C
	rm -f main.bbl main.run.xml
