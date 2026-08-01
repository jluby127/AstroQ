# Build the HIRES 2026B paper.
#
#   make          build main.pdf
#   make params   regenerate params.tex from data/
#   make check    fail if any \param{} key is unresolved in the PDF
#   make clean    remove latexmk build products

PYTHON ?= python
LATEXMK ?= latexmk

# Run directory holding outputs/semester_planner.h5 for the first night of the
# semester. Only needed by `make figures`; the PDFs it produces are committed
# so the paper builds without it.
RUN_DIR ?= ../../2026B/2026-08-01/band1_baseline

.PHONY: all params tables figures check clean

all: main.pdf

# params.tex and the tables are generated, never hand-edited, so they are
# rebuilt whenever any recorded artifact changes.
params: params.tex
tables: table_semester.tex table_programs.tex

params.tex: tools/make_params.py $(wildcard data/*/*)
	$(PYTHON) tools/make_params.py

table_semester.tex table_programs.tex: tools/make_tables.py data/01/baseline-report-final.txt
	$(PYTHON) tools/make_tables.py

# Figures need a solved run directory, so they are not wired into main.pdf.
# Regenerate deliberately when the underlying plan changes.
figures:
	cd ../.. && PYTHONPATH=. $(PYTHON) papers/hires-2026b/tools/make_figures.py $(RUN_DIR:../../%=%)

main.pdf: main.tex params.tex table_semester.tex table_programs.tex variablesTable.tex bib.bib
	$(LATEXMK) -pdf -interaction=nonstopmode -halt-on-error main.tex

# An unresolved \param{} key renders as a literal red "XX" in the PDF. Catch it
# at the source level instead, which needs no pdftotext and names the key.
check:
	$(PYTHON) tools/make_params.py --check

clean:
	$(LATEXMK) -C
	rm -f main.bbl main.run.xml
