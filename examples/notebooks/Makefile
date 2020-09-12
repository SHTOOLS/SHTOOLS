###############################################################################
#
#   make all
#       Convert notebooks to html.
#
#   make clean
#       Remove html files.
#
#   make run-notebooks
#       Run notebooks and check for errors
#
###############################################################################

JUPYTER ?= jupyter nbconvert
PYTHON = python3

SHELL=/bin/sh

WWWNBDIR = ../../docs/pages/mydoc/notebooks
JFLAGS = --to html --execute --ExecutePreprocessor.timeout=120 --output-dir=$(WWWNBDIR)

IPYNB = grids-and-coefficients.ipynb \
	localized-spectral-analysis.ipynb \
	gravity-and-magnetic-fields.ipynb \
	plotting-maps.ipynb \
	low-level-spherical-harmonic-analyses.ipynb \
	advanced-localized-spectral-analysis.ipynb \
	advanced-shcoeffs-and-shgrid-usage.ipynb \
	spherical-harmonic-normalizations.ipynb \
	advanced-shwindow-usage.ipynb \
	3d-plots.ipynb

HTML = $(addprefix $(WWWNBDIR)/, grids-and-coefficients.html \
	localized-spectral-analysis.html \
	gravity-and-magnetic-fields.html \
	plotting-maps.html \
	low-level-spherical-harmonic-analyses.html \
	advanced-localized-spectral-analysis.html \
	advanced-shcoeffs-and-shgrid-usage.html \
	spherical-harmonic-normalizations.html \
	advanced-shwindow-usage.html \
	3d-plots.html)

.PHONY: clean all run-notebooks

all: $(HTML)

clean:
	@-rm -f $(HTML)

run-notebooks:
	@$(PYTHON) test_notebooks.py

$(WWWNBDIR)/%.html: %.ipynb
	$(JUPYTER) $(JFLAGS) $<
