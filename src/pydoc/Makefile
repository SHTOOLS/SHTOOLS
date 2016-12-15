###############################################################################
#
#   make all
#       Make man and html-man pages, and then move to their appropriate
#       locations.
#
#   make install-man
#       Make man pages from markdown source files.
#
#   make install-html
#       Make html-man pages from markdown source files.
#
#   make clean
#       Remove the man and html-man pages in the LOCAL directory.
#
#   make uninstall
#       Remove man pages from /usr/local
#
###############################################################################

SHELL=/bin/sh

PANDOC = pandoc
PDHTMLFLAGS = --base-header-level=1
MANDIR = ../../man/man1
WWWDIR = ../../www/man/python

MDFILES = pyplmindex.md pyplmbar.md pyplmbar_d1.md pyplbar.md pyplbar_d1.md \
	pyplon.md pyplmon.md pyplmon_d1.md pyplon_d1.md pyplmschmidt.md \
	pyplschmidt.md pyplmschmidt_d1.md pyplschmidt_d1.md pyplegendrea.md \
	pyplegendre.md pyplegendrea_d1.md pyplegendre_d1.md pyshexpanddh.md \
	pymakegriddh.md pyshexpanddhc.md pymakegriddhc.md pyshglq.md \
	pyglqgridcoord.md pyshexpandglq.md pymakegridglq.md pyshexpandglqc.md \
	pymakegridglqc.md pyshexpandlsq.md pymakegrid2d.md pymakegridpoint.md \
	pymakegridpointc.md pyshmultiply.md pyshread.md pyshreaderror.md \
	pyshreadh.md pyshreaderrorh.md pyshread2.md pyshread2error.md \
	pyshreadjpl.md pyshreadjplerror.md pyshcilmtocindex.md \
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
	pynormalgravity.md pymakemaggriddh.md pyshmagpowerspectrum.md \
	pyshmagpowerl.md pyshmtcouplingmatrix.md pyshmultitapermaskcse.md \
	pyshmultitapermaskse.md pyshbiaskmask.md read_icgem_gfc.md \
	spectrum.md cross_spectrum.md pydhaj.md

MANFILES = $(addprefix $(MANDIR)/, pyplmindex.1 pyplmbar.1 pyplmbar_d1.1 \
	pyplbar.1 pyplbar_d1.1 pyplon.1 pyplmon.1 pyplmon_d1.1 pyplon_d1.1 \
	pyplmschmidt.1 pyplschmidt.1 pyplmschmidt_d1.1 pyplschmidt_d1.1 \
	pyplegendrea.1 pyplegendre.1 pyplegendrea_d1.1 pyplegendre_d1.1 \
	pyshexpanddh.1 pymakegriddh.1 pyshexpanddhc.1 pymakegriddhc.1 pyshglq.1 \
	pyglqgridcoord.1 pyshexpandglq.1 pymakegridglq.1 pyshexpandglqc.1 \
	pymakegridglqc.1 pyshexpandlsq.1 pymakegrid2d.1 pymakegridpoint.1 \
	pymakegridpointc.1 pyshmultiply.1 pyshread.1 pyshreaderror.1 \
	pyshreadh.1 pyshreaderrorh.1 pyshread2.1 pyshread2error.1 pyshreadjpl.1 \
	pyshreadjplerror.1 pyshcilmtocindex.1 pyshcindextocilm.1 \
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
	pymakemaggriddh.1 pyshmagpowerspectrum.1 pyshmagpowerl.1 \
	pyshmtcouplingmatrix.1 pyshmultitapermaskcse.1 pyshmultitapermaskse.1 \
	pyshbiaskmask.1 read_icgem_gfc.1 spectrum.1 cross_spectrum.1 \
	pydhaj.1)

HTMLFILES_legendre = $(addprefix $(WWWDIR)/, pyplmindex.html pyplmbar.html \
	pyplmbar_d1.html pyplbar.html pyplbar_d1.html pyplon.html pyplmon.html \
	pyplmon_d1.html pyplon_d1.html pyplmschmidt.html pyplschmidt.html \
	pyplmschmidt_d1.html pyplschmidt_d1.html pyplegendrea.html \
	pyplegendre.html pyplegendrea_d1.html pyplegendre_d1.html)

HTMLFILES_transformations = $(addprefix $(WWWDIR)/, pyshexpanddh.html \
	pymakegriddh.html pyshexpanddhc.html pymakegriddhc.html pyshglq.html \
	pyglqgridcoord.html pyshexpandglq.html pymakegridglq.html \
	pyshexpandglqc.html pymakegridglqc.html pyshexpandlsq.html \
	pymakegrid2d.html pymakegridpoint.html pymakegridpointc.html \
	pyshmultiply.html)

HTMLFILES_io = $(addprefix $(WWWDIR)/, pyshread.html pyshreaderror.html \
	pyshreadh.html pyshreaderrorh.html pyshread2.html pyshread2error.html \
	pyshreadjpl.html pyshreadjplerror.html pyshcilmtocindex.html \
	pyshcindextocilm.html pyshcilmtovector.html pyshvectortocilm.html \
	pyshrtoc.html pyshctor.html pyyilmindexvector.html read_icgem_gfc.html)

HTMLFILES_spectralanalysis = $(addprefix $(WWWDIR)/, spectrum.html \
	cross_spectrum.html pyshadmitcorr.html pyshconfidence.html \
	pyshmultitaperse.html \
	pyshmultitapercse.html pyshlocalizedadmitcorr.html pyshreturntapers.html \
	pyshreturntapersm.html pycomputedm.html pycomputedg82.html \
	pyshfindlwin.html pyshbiask.html pyshbiasadmitcorr.html pyshmtdebias.html \
	pyshmtvaropt.html pyshsjkpg.html pyshreturntapersmap.html \
	pycomputedmap.html pycurve2mask.html pyshbias.html \
	pysphericalcapcoef.html pyshmtcouplingmatrix.html \
	pyshmultitapermaskcse.html pyshmultitapermaskse.html pyshbiaskmask.html)

HTMLFILES_rotations = $(addprefix $(WWWDIR)/, pydjpi2.html \
	pyshrotatecoef.html pyshrotaterealcoef.html)

HTMLFILES_gravmag = $(addprefix $(WWWDIR)/, pymakegravgriddh.html \
	pymakegravgradgriddh.html pymakegeoidgriddh.html pycilmplusdh.html \
	pycilmminusdh.html pycilmplusrhohdh.html pycilmminusrhohdh.html \
	pybatohilmdh.html pybatohilmrhohdh.html pydowncontfilterma.html \
	pydowncontfiltermc.html pynormalgravity.html pymakemaggriddh.html \
	pyshmagpowerl.html pyshmagpowerspectrum.html)

HTMLFILES_utils = $(addprefix $(WWWDIR)/, pymakecirclecoord.html \
	pymakeellipsecoord.html pywigner3j.html pydhaj.html)


all: install-man install-html


install-man: $(MANFILES)
	@echo
	@echo "--> Python man pages created successfully"
	@echo

install-html: $(HTMLFILES_modules) $(HTMLFILES_legendre) $(HTMLFILES_transformations) \
	$(HTMLFILES_io) $(HTMLFILES_spectralanalysis) \
	$(HTMLFILES_rotations) $(HTMLFILES_gravmag) \
	$(HTMLFILES_utils)
	@echo
	@echo "--> Python html web pages created successfully"
	@echo

$(HTMLFILES_legendre) : $(WWWDIR)/%.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_legendre.html test.html tail_legendre.html > $@
	@rm -f test.html

$(HTMLFILES_transformations) : $(WWWDIR)/%.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_transformations.html test.html tail_transformations.html > $@
	@rm -f test.html

$(HTMLFILES_io) : $(WWWDIR)/%.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_io.html test.html tail_io.html > $@
	@rm -f test.html

$(HTMLFILES_spectralanalysis) : $(WWWDIR)/%.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_spectralanalysis.html test.html tail_spectralanalysis.html > $@
	@rm -f test.html

$(HTMLFILES_rotations) : $(WWWDIR)/%.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_rotations.html test.html tail_rotations.html > $@
	@rm -f test.html

$(HTMLFILES_gravmag) : $(WWWDIR)/%.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_gravmag.html test.html tail_gravmag.html > $@
	@rm -f test.html

$(HTMLFILES_utils) : $(WWWDIR)/%.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_utils.html test.html tail_utils.html > $@
	@rm -f test.html


.PHONY: clean uninstall

clean:
	-rm -f $(MANFILES) $(HTMLFILES) 

uninstall:
	-cd /usr/local/share/man/man1/ ; rm $(MANFILES)


.SUFFIXES: $(SUFFIXES) .md .html .1


$(MANDIR)/%.1 : %.md
	printf "%% $(basename $<)(1) Python | SHTOOLS $(VERSION)\n%%\n%% %s" > temp.md
	stat -t "%F" -f "%Sm" $< >> temp.md
	cat $< >> temp.md
	$(PANDOC) -s -t man temp.md -o $@ --filter delink.hs
	rm -f temp.md

