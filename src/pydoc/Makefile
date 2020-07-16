###############################################################################
#
#   make all
#       First, run the python script doc2md that converts the pure python
#       function docstrings into md files. Then, run make_www_doc.py that
#       reformats the md documentation files for use with the web documentation
#       found in doc/pages/mydoc/pydoc.
#
###############################################################################

SHELL=/bin/sh

PYTHON = python3

all:
	$(PYTHON) doc2md.py
	$(PYTHON) make_www_doc.py
