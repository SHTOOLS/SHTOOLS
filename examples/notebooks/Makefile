###############################################################################
#
#   make all
#       Run notebooks and convert to html.
#
#   make clean
#       Remove html files.
#
###############################################################################

JUPYTER ?= jupyter nbconvert

SHELL=/bin/sh

WWWNBDIR = ../../docs/pages/mydoc/notebooks
JFLAGS = --to html --execute --ExecutePreprocessor.timeout=120 --output-dir=$(WWWNBDIR)

IPYNB = grids-and-coefficients.ipynb \
	localized-spectral-analysis.ipynb \
	gravity-and-magnetic-fields.ipynb \
	plotting-maps.ipynb
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

all: $(HTML)

.PHONY: clean all

clean:
	@-rm -f $(HTML)

$(WWWNBDIR)/%.html: %.ipynb
	$(JUPYTER) $(JFLAGS) $<
