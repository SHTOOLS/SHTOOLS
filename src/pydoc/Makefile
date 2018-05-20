###############################################################################
#
#   make all
#       Make man pages, and then move to their appropriate locations.
#
#   make install-man
#       Make man pages from markdown source files.
#
#   make clean
#       Remove the man pages in the LOCAL directory.
#
#   make uninstall
#       Remove man pages from /usr/local
#
###############################################################################

SHELL=/bin/sh

PANDOC = pandoc
MANDIR = ../../man/man1
PYTHON = python3

MDFILES = pyplmindex.md pyplmbar.md pyplmbar_d1.md pyplbar.md pyplbar_d1.md \
	pyplon.md pyplmon.md pyplmon_d1.md pyplon_d1.md pyplmschmidt.md \
	pyplschmidt.md pyplmschmidt_d1.md pyplschmidt_d1.md pyplegendrea.md \
	pyplegendre.md pyplegendrea_d1.md pyplegendre_d1.md pyshexpanddh.md \
	pymakegriddh.md pyshexpanddhc.md pymakegriddhc.md pyshglq.md \
	pyglqgridcoord.md pyshexpandglq.md pymakegridglq.md pyshexpandglqc.md \
	pymakegridglqc.md pyshexpandlsq.md pymakegrid2d.md pymakegridpoint.md \
	pymakegridpointc.md pyshmultiply.md pyshread.md pyshread2.md \
	pyshread2error.md pyshreadjpl.md pyshreadjplerror.md pyshcilmtocindex.md \
	pyshcindextocilm.md pyshcilmtovector.md pyshvectortocilm.md pyshrtoc.md \
	pyshctor.md pyshadmitcorr.md pyshconfidence.md pyshmultitaperse.md \
	pyshmultitapercse.md pyshlocalizedadmitcorr.md pyshreturntapers.md \
	pyshreturntapersm.md pycomputedm.md pycomputedg82.md pyshfindlwin.md \
	pyshbiask.md pyshbiasadmitcorr.md pyshmtdebias.md pyyilmindexvector.md \
	pyshmtvaropt.md pyshsjkpg.md pyshreturntapersmap.md pycomputedmap.md \
	pycurve2mask.md pyshbias.md pysphericalcapcoef.md pydjpi2.md \
	pyshrotatecoef.md pyshrotaterealcoef.md pymakecirclecoord.md \
	pymakeellipsecoord.md pywigner3j.md pymakegravgriddh.md \
	pymakegravgradgriddh.md pymakegeoidgriddh.md pycilmplusdh.md \
	pycilmminusdh.md pycilmplusrhohdh.md pycilmminusrhohdh.md pybatohilmdh.md \
	pybatohilmrhohdh.md pydowncontfilterma.md pydowndontfiltermc.md \
	pynormalgravity.md pymakemaggriddh.md pyshmtcouplingmatrix.md \
	pyshmultitapermaskcse.md pyshmultitapermaskse.md pyshbiaskmask.md \
	read_icgem_gfc.md spectrum.md cross_spectrum.md pydhaj.md mag_spectrum.md \
	convert.md

MANFILES = $(addprefix $(MANDIR)/, pyplmindex.1 pyplmbar.1 pyplmbar_d1.1 \
	pyplbar.1 pyplbar_d1.1 pyplon.1 pyplmon.1 pyplmon_d1.1 pyplon_d1.1 \
	pyplmschmidt.1 pyplschmidt.1 pyplmschmidt_d1.1 pyplschmidt_d1.1 \
	pyplegendrea.1 pyplegendre.1 pyplegendrea_d1.1 pyplegendre_d1.1 \
	pyshexpanddh.1 pymakegriddh.1 pyshexpanddhc.1 pymakegriddhc.1 pyshglq.1 \
	pyglqgridcoord.1 pyshexpandglq.1 pymakegridglq.1 pyshexpandglqc.1 \
	pymakegridglqc.1 pyshexpandlsq.1 pymakegrid2d.1 pymakegridpoint.1 \
	pymakegridpointc.1 pyshmultiply.1 pyshread.1 pyshread2.1 pyshread2error.1 \
	pyshreadjpl.1 pyshreadjplerror.1 pyshcilmtocindex.1 pyshcindextocilm.1 \
	pyshcilmtovector.1 pyshvectortocilm.1 pyshrtoc.1 pyshctor.1 \
	pyshadmitcorr.1 pyshconfidence.1 pyshmultitaperse.1 pyshmultitapercse.1 \
	pyshlocalizedadmitcorr.1 pyshreturntapers.1 pyshreturntapersm.1 \
	pycomputedm.1 pycomputedg82.1 pyshfindlwin.1 pyshbiask.1 \
	pyshbiasadmitcorr.1 pyshmtdebias.1 pyyilmindexvector.1 pyshmtvaropt.1 \
	pyshsjkpg.1 pyshreturntapersmap.1 pycomputedmap.1 pycurve2mask.1 \
	pyshbias.1 pysphericalcapcoef.1 pydjpi2.1 pyshrotatecoef.1 \
	pyshrotaterealcoef.1 pymakecirclecoord.1 pymakeellipsecoord.1\
	pywigner3j.1 pymakegravgriddh.1 pymakegravgradgriddh.1 \
	pymakegeoidgriddh.1 pycilmplusdh.1 pycilmminusdh.1 pycilmplusrhohdh.1 \
	pycilmminusrhohdh.1 pybatohilmdh.1 pybatohilmrhohdh.1 \
	pydowncontfilterma.1 pydowncontfiltermc.1 pynormalgravity.1 \
	pymakemaggriddh.1 pyshmtcouplingmatrix.1 pyshmultitapermaskcse.1 \
	pyshmultitapermaskse.1 pyshbiaskmask.1 read_icgem_gfc.1 spectrum.1 \
	cross_spectrum.1 pydhaj.1 mag_spectrum.1 convert.1)


all: install-man
	$(PYTHON) make_www_doc.py

install-man: $(MANFILES)
	@echo
	@echo "--> Python man pages created successfully"
	@echo


.PHONY: clean uninstall

clean:
	-rm -f $(MANFILES)

uninstall:
	-cd /usr/local/share/man/man1/ ; rm $(MANFILES)


.SUFFIXES: $(SUFFIXES) .md .1


$(MANDIR)/%.1 : %.md
	printf "%% $(basename $<)(1) Python | SHTOOLS $(VERSION)\n%%\n%% %s" > temp.md
	stat -t "%F" -f "%Sm" $< >> temp.md
	cat $< >> temp.md
	$(PANDOC) -s -t man temp.md -o $@ --filter delink.hs
	rm -f temp.md

