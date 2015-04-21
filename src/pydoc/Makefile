######################################################################
#
#	make all
#		Make man and html-man pages, and then move to their
#		appropriate locations.
#
#	make install-man
#		Make man pages from markdown source files.
#	
#	make install-html
#		Make html-man pages from markdown source files.
#
#	make clean
#		Remove the man and html-man pages in the LOCAL directory.
#
#	make uninstall
#		Remove man pages from /usr/local
#
#	Written by Mark Wieczorek (July 2006).
#
#######################################################################

SHELL=/bin/tcsh

PANDOC = pandoc
PDHTMLFLAGS = --base-header-level=1
MANDIR = ../../man/man1
WWWDIR = ../../www/man/python

MDFILES = pyplmindex.md pyplmbar.md pyplmbar_d1.md pyplbar.md pyplbar_d1.md pyplon.md \
    pyplmon.md pyplmon_d1.md pyplon_d1.md pyplmschmidt.mdpyplschmidt.md \
    pyplmschmidt_d1.md pyplschmidt_d1.md pyplegendrea.md pyplegendre.md pyplegendrea_d1.md \
    pyplegendre_d1.md pyshexpanddh.md pymakegriddh.md pyshexpanddhc.md pymakegriddhc.md \
    pyshglq.md pyglqgridcoord.md pyshexpandglq.md pymakegridglq.md pyshexpandglqc.md \
    pymakegridglqc.md pyshexpandlsq.md pymakegrid2d.md pymakegridpoint.md pyshmultiply.md \
    pyshread.md pyshreaderror.md pyshreadh.md pyshreaderrorh.md pyshread2.md \
    pyshread2error.md pyshreadjpl.md pyshreadjplerror.md pyshcilmtocindex.md \
    pyshcindextocilm.md pyshcilmtovector.md pyshvectortocilm.md pyshrtoc.md \
    pyshctor.md pyshpowerl.md pyshpowerdensityl.md pyshcrosspowerl.md \
    pyshcrosspowerdensityl.md pyshpowerspectrum.md pyshpowerspectrumdensity.md \
    pyshcrosspowerspectrum.md pyshcrosspowerspectrumdensity.md pyshadmitcorr.md \
    pyshconfidence.md pyshpowerlc.md pyshpowerdensitylc.md pyshcrosspowerlc.md \
    pyshcrosspowerdensitylc.md pyshpowerspectrumc.md pyshpowerspectrumdensityc.md \
    pyshcrosspowerspectrumc.md pyshcrosspowerspectrumdensityc.md pyshmultitaperse.md \
    pyshmultitapercse.md pyshlocalizedadmitcorr.md pyshreturntapers.md pyshreturntapersm.md \
    pycomputedm.md pycomputedg82.md pyshfindlwin.md pyshbiask.md pyshbiasadmitcorr.md \
    pyshmtdebias.md pyyilmindexvector.md
#pycilmplusdh.pod pycilmplusglq.pod pycilmminusdh.pod pycilmminusglq.pod   \
	 pydjpi2.pod pyeigvalsym.pod pyeigvalvecsym.pod pyeigvalvecsymtri.pod \
	pyhilmdh.pod pyhilmglq.pod pymakecirclecoord.pod pymakegeoidgriddh.pod pymakegeoidgridglq.pod \
	pymakegeoidgrid2d.pod pymakemaggriddh.pod pypreglq.pod \
	 pyrandomn.pod pyrandomgaussian.pod pyshbias.pod     \
	 pyshmtvaropt.pod   \
	pywl.pod pywl_curv.pod    \
	   pyshrotatecoef.pod pyshrotaterealcoef.pod \
	pywigner3j.pod pyshsjkpg.pod  pysphericalcapcoef.pod  \
	pymakegravgriddh.pod  pyshmagpowerl.pod pyshmagpowerspectrum.pod     \
	 pycomputedmap.pod pycurve2mask.pod pymakeellipsecoord.pod \
	pycilmplusrhohdh.pod pycilmminusrhohdh.pod pycilmplusrhohglq.pod pycilmminusrhohglq.pod pyhilmrhohdh.pod pyhilmrhohglq.pod pynormalgravity.pod \
	pymakegravgradgriddh.pod 

MANFILES = $(addprefix $(MANDIR)/, pyplmindex.1 pyplmbar.1 pyplmbar_d1.1 pyplbar.1 pyplbar_d1.1 \
    pyplon.1 pyplmon.1 pyplmon_d1.1 pyplon_d1.1 pyplmschmidt.1 pyplschmidt.1 pyplmschmidt_d1.1 \
    pyplschmidt_d1.1 pyplegendrea.1 pyplegendre.1 pyplegendrea_d1.1 pyplegendre_d1.1 \
    pyshexpanddh.1 pymakegriddh.1 pyshexpanddhc.1 pymakegriddhc.1 pyshglq.1 pyglqgridcoord.1\
    pyshexpandglq.1 pymakegridglq.1 pyshexpandglqc.1 pymakegridglqc.1 pyshexpandlsq.1 \
    pymakegrid2d.1 pymakegridpoint.1 pyshmultiply.1 pyshread.1 pyshreaderror.1 \
    pyshreadh.1 pyshreaderrorh.1 pyshread2.1 pyshread2error.1 pyshreadjpl.1 \
    pyshreadjplerror.1 pyshcilmtocindex.1 pyshcindextocilm.1 pyshcilmtovector.1 \
    pyshvectortocilm.1 pyshrtoc.1 pyshctor.1 pyshpowerl.1 pyshpowerdensityl.1 \
    pyshcrosspowerl.1 pyshcrosspowerdensityl.1 pyshpowerspectrum.1 \
    pyshpowerspectrumdensity.1 pyshcrosspowerspectrum.1 pyshcrosspowerspectrumdensity.1 \
    pyshadmitcorr.1 pyshconfidence.1 pyshpowerlc.1 pyshpowerdensitylc.1 \
    pyshcrosspowerlc.1 pyshcrosspowerdensitylc.1 pyshpowerspectrumc.1 \
    pyshpowerspectrumdensityc.1 pyshcrosspowerspectrumc.1 pyshcrosspowerspectrumdensityc.1 \
    pyshmultitaperse.1 pyshmultitapercse.1 pyshlocalizedadmitcorr.1 pyshreturntapers.1 \
    pyshreturntapersm.1 pycomputedm.1 pycomputedg82.1 pyshfindlwin.1 pyshbiask.1 \
    pyshbiasadmitcorr.1 pyshmtdebias.1 pyyilmindexvector.1)
#pycilmplusdh.1 pycilmplusglq.1 pycilmminusdh.1 pycilmminusglq.1   \
	 pydjpi2.1 pyeigvalsym.1 pyeigvalvecsym.1 pyeigvalvecsymtri.1 pyhilmdh.1 pyhilmglq.1 \
	pymakecirclecoord.1 pymakegeoidgriddh.1 pymakegeoidgridglq.1 pymakegeoidgrid2d.1 \
	  pymakemaggriddh.1  pypreglq.1 pyrandomn.1 pyrandomgaussian.1 \
	 pyshbias.1     \
	 pyshmtvaropt.1   pywl.1 pywlcurv.1     \
	  pyshrotatecoef.1 pyshrotaterealcoef.1 \
	pywigner3j.1 pyshsjkpg.1  pysphericalcapcoef.1  pymakegravgriddh.1 \
	 pyshmagpowerl.1 pyshmagpowerspectrum.1      \
	pycomputedmap.1 pyshreturntapersmap.1 pycurve2mask.1 pymakeellipsecoord.1 pycilmplusrhohdh.1  pycilmminusrhohdh.1 pycilmplusrhohglq.1 pycilmminusrhohglq.1 \
	pyhilmrhohdh.1 pyhilmrhohglq.1 pynormalgravity.1 pymakegravgradgriddh.1 

HTMLFILES_legendre =  $(addprefix $(WWWDIR)/, pyplmindex.html pyplmbar.html pyplmbar_d1.html pyplbar.html \
    pyplbar_d1.html pyplon.html pyplmon.html pyplmon_d1.html pyplon_d1.html pyplmschmidt.html \
    pyplschmidt.html pyplmschmidt_d1.html pyplschmidt_d1.html pyplegendrea.html pyplegendre.html \
    pyplegendrea_d1.html pyplegendre_d1.html )
	  
HTMLFILES_transformations = $(addprefix $(WWWDIR)/, pyshexpanddh.html pymakegriddh.html  \
	pyshexpanddhc.html pymakegriddhc.html pyshglq.html pyglqgridcoord.html pyshexpandglq.html\
    pymakegridglq.html pyshexpandglqc.html pymakegridglqc.html pyshexpandlsq.html \
    pymakegrid2d.html pymakegridpoint.html pyshmultiply.html)

HTMLFILES_io = $(addprefix $(WWWDIR)/, pyshread.html pyshreaderror.html pyshreadh.html \
    pyshreaderrorh.html pyshread2.html pyshread2error.html pyshreadjpl.html \
    pyshreadjplerror.html pyshcilmtocindex.html pyshcindextocilm.html pyshcilmtovector.html \
    pyshvectortocilm.html pyshrtoc.html pyshctor.html pyyilmindexvector.html)
 	  
HTMLFILES_gspectra = $(addprefix $(WWWDIR)/, pyshpowerl.html pyshpowerdensityl.html \
	pyshcrosspowerl.html pyshcrosspowerdensityl.html pyshpowerspectrum.html \
	pyshpowerspectrumdensity.html pyshcrosspowerspectrum.html \
	pyshcrosspowerspectrumdensity.html pyshadmitcorr.html pyshconfidence.html \
	pyshpowerlc.html pyshpowerdensitylc.html pyshcrosspowerlc.html \
	pyshcrosspowerdensitylc.html pyshpowerspectrumc.html pyshpowerspectrumdensityc.html \
	pyshcrosspowerspectrumc.html pyshcrosspowerspectrumdensityc.html)

HTMLFILES_localized = $(addprefix $(WWWDIR)/, pyshmultitaperse.html \
	pyshmultitapercse.html pyshlocalizedadmitcorr.html pyshreturntapers.html \
	pyshreturntapersm.html pycomputedm.html pycomputedg82.html pyshfindlwin.html \
	pyshbiask.html pyshbiasadmitcorr.html pyshmtdebias.html)
#  pyshbias.html   \
	 pyshmtvaropt.html   \
	  pyshsjkpg.html \
	 pysphericalcapcoef.html   \
	pycomputedmap.html pyshreturntapersmap.html pycurve2mask.html

HTMLFILES_rotations = 
#pydjpi2.html pyshrotatecoef.html pyshrotaterealcoef.html

HTMLFILES_gravmag = 
#pycilmplusdh.html pycilmplusglq.html pycilmminusdh.html pycilmminusglq.html pyhilmdh.html pyhilmglq.html pymakegeoidgriddh.html \
	pymakegeoidgridglq.html pymakegeoidgrid2d.html pymakemaggriddh.html pywl.html pywlcurv.html pymakegravgriddh.html \
	pyshmagpowerl.html pyshmagpowerspectrum.html pycilmplusrhohdh.html pycilmplusrhohglq.html pycilmminusrhohdh.html pycilmminusrhohglq.html pyhilmrhohdh.html \
	pyhilmrhohglq.html pynormalgravity.html pymakegravgradgriddh.html

HTMLFILES_other =  
#pyeigvalsym.html pyeigvalvecsym.html pyeigvalvecsymtri.html \
	pymakecirclecoord.html pyrandomn.html pyrandomgaussian.html \
	pywigner3j.html pymakeellipsecoord.html pypreglq.html
 

all: install-man install-html


install-man: $(MANFILES)
	@echo
	@echo PYTHON MAN PAGES SUCCESSFULLY MADE
	@echo

install-html: $(HTMLFILES_modules) $(HTMLFILES_legendre) $(HTMLFILES_transformations) \
	$(HTMLFILES_io) $(HTMLFILES_gspectra) $(HTMLFILES_localized) \
	$(HTMLFILES_rotations) $(HTMLFILES_gravmag) \
	$(HTMLFILES_other)
	@echo
	@echo PYTHON HTML WEB PAGES SUCCESSFULLY MADE
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

$(HTMLFILES_gspectra) : $(WWWDIR)/%.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_gspectra.html test.html tail_gspectra.html > $@
	@rm -f test.html

$(HTMLFILES_localized) : $(WWWDIR)/%.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_localized.html test.html tail_localized.html > $@
	@rm -f test.html

$(HTMLFILES_rotations) : $(WWWDIR)/%.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_rotations.html test.html tail_rotations.html > $@
	@rm -f test.html

$(HTMLFILES_gravmag) : $(WWWDIR)/%.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_gravmag.html test.html tail_gravmag.html > $@
	@rm -f test.html

$(HTMLFILES_other) : $(WWWDIR)/%.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_other.html test.html tail_other.html > $@
	@rm -f test.html


.PHONY: clean uninstall

clean:
	-rm -f $(MANFILES) $(HTMLFILES) 

uninstall:
	-cd /usr/local/share/man/man1/ ; rm $(MANFILES)


.SUFFIXES: $(SUFFIXES)  .md .html .1
	

$(MANDIR)/%.1 : %.md
	printf "%% $(basename $<)(1) SHTOOLS $(VERSION) | SHTOOLS $(VERSION)\n%%\n\%% %s" > temp.md
	stat -t "%F" -f "%Sm" $< >> temp.md
	cat $< >> temp.md
	$(PANDOC) -s -t man temp.md -o $@ --filter delink.hs
	rm -f temp.md

