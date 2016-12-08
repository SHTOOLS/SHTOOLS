###############################################################################
#
#   make all
#       Run notebooks and convert to html.
#
#   make clean
#       Remove html files.
#
###############################################################################

JUPYTER = jupyter nbconvert

SHELL=/bin/sh

WWWNBDIR = ../../www/notebooks
JFLAGS = --to html --execute --ExecutePreprocessor.timeout=120 --output-dir=$(WWWNBDIR)

IPYNB = Introduction-1.ipynb Introduction-2.ipynb tutorial_1.ipynb \
	tutorial_2.ipynb tutorial_3.ipynb tutorial_4.ipynb tutorial_5.ipynb \
	tutorial_6.ipynb

HTML = $(addprefix $(WWWNBDIR)/, Introduction-1.html Introduction-2.html \
	tutorial_1.html tutorial_2.html tutorial_3.html tutorial_4.html \
	tutorial_5.html tutorial_6.html)

all: $(HTML)

.PHONY: clean all

clean:
	@-rm -f $(HTML)

$(WWWNBDIR)/%.html: %.ipynb
	$(JUPYTER) $(JFLAGS) $<
