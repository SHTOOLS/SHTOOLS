######################################################################
#
#	make all
#		Make man and html-man pages, and then move to their
#		appropriate locations.
#
#	make install-man
#		Make man pages from POD files.
#	
#	make install-html
#		Make html-man pages from the previously created man pages.
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
PDHTMLFLAGS = --base-header-level=2
MANDIR = ../../man/man1
WWWDIR = ../../www/man/python

MDFILES = pymakegriddh.md pyplmbar.pod
#pycilmplusdh.pod pycilmplusglq.pod pycomputed0.pod pycomputedg82.pod pycomputedm.pod \
	pydhaj.pod pydjpi2.pod pyeigvalsym.pod pyeigvalvecsym.pod pyeigvalvecsymtri.pod \
	pyhilmdh.pod pyhilmglq.pod pymakecirclecoord.pod pymakegeoidgriddh.pod pymakegeoidgridglq.pod \
	pymakegeoidgrid2d.pod pymakegrid2d.pod pyglqgridcoord.pod \
	pymakegridglq.pod pymakegridpoint.pod pymakemaggriddh.pod \
 	pyplbar_d1.pod pyplbar.pod pyplegendre.pod pyplegendre_d1.pod \
	pyplegendrea_d1.pod pyplegendrea.pod pyplmbar_d1.pod  pyplmindex.pod \
	pyplon.pod pyplon_d1.pod pyplmon.pod pyplmon_d1.pod pyplschmidt.pod pyplschmidt_d1.pod \
	pyplmschmidt.pod pyplmschmidt_d1.pod pyprecompute.pod pypreglq.pod pynglq.pod pynglqsh.pod \
	pynglqshn.pod pyrandomn.pod pyrandomgaussian.pod \
	pyshadmitcorr.pod pyshbias.pod pyshbiask.pod pyshrtoc.pod \
	pyshctor.pod pyshcilmtocindex.pod pyshcindextocilm.pod  \
	pyshexpanddh.pod pyshexpandglq.pod pyshexpandlsq.pod pyshfindlwin.pod \
	pyshlocalizedadmitcorr.pod pyshmtvaropt.pod pyshmtvaropt0.pod pyshmultiply.pod \
	pywl.pod pywl_curv.pod pyshmultitapercse.pod pyshmultitaperse.pod pyshpowerl.pod \
	pyshpowerdensityl.pod pyshcrosspowerdensityl.pod pyshpowerspectrum.pod \
	pyshpowerspectrumdensity.pod pyshcrosspowerspectrum.pod pyshcrosspowerspectrumdensity.pod \
	pyshcrosspowerl.pod pyshread.pod pyshreaderror.pod pyshreadjpl.pod pyshreadjplerror.pod pyshread2.pod \
	pyshread2error.pod pyshreturntapers.pod pyshreturntapersm.pod pyshrotatecoef.pod pyshrotaterealcoef.pod \
	pywigner3j.pod pyshsjkpg.pod pyshsjkpg0.pod pysphericalcapcoef.pod pyshmtdebias.pod \
	pymakegravgriddh.pod pyshconfidence.pod pyshmagpowerl.pod pyshmagpowerspectrum.pod \
	pymakegriddhc.pod pyshexpanddhc.pod pymakegridglqc.pod pyshexpandglqc.pod pyshpowerlc.pod \
	pyshpowerdensitylc.pod pyshcrosspowerdensitylc.pod pyshpowerspectrumc.pod \
	pyshpowerspectrumdensityc.pod pyshcrosspowerspectrumc.pod pyshcrosspowerspectrumdensityc.pod \
	pyshcrosspowerlc.pod pyshbiasadmitcorr.pod pyshcilmtovector.pod pyshvectortocilm.pod \
	pyyilmindex.pod pycomputedmap.pod pyshreturntapers.pod pycurve2mask.pod pymakeellipsecoord.pod \
	pycilmplusrhohdh.pod pycilmplusglq.pod pyhilmrhohdh.pod pyhilmrhohglq.pod pynormalgravity.pod \
	pymakegravgradgriddh.pod pyshreadh.pod pyshreaderrorh.pod

MANFILES = pymakegriddh.1 pyplmbar.1
#pycilmplusdh.1 pycilmplusglq.1 pycomputed0.1 pycomputedg82.1 pycomputedm.1 \
	pydhaj.1 pydjpi2.1 pyeigvalsym.1 pyeigvalvecsym.1 pyeigvalvecsymtri.1 pyhilmdh.1 pyhilmglq.1 \
	pymakecirclecoord.1 pymakegeoidgriddh.1 pymakegeoidgridglq.1 pymakegeoidgrid2d.1 \
	pymakegrid2d.1 pyglqgridcoord.1  \
	pymakegridglq.1 pymakegridpoint.1 pymakemaggriddh.1 pyplbar_d1.1 \
	pyplbar.1 pyplegendre.1 pyplegendre_d1.1 pyplegendrea_d1.1 pyplegendrea.1 pyplmbar_d1.1 \
	 pyplmindex.1 pyplon.1 pyplon_d1.1 pyplmon.1 pyplmon_d1.1 pyplschmidt.1 pyplschmidt_d1.1 \
	pyplmschmidt.1 pyplmschmidt_d1.1 pyprecompute.1 pypreglq.1 pynglq.1 pynglqsh.1 pynglqshn.1 \
	pyrandomn.1 pyrandomgaussian.1 \
	pyshadmitcorr.1 pyshbias.1 pyshbiask.1 pyshrtoc.1 pyshctor.1 pyshcilmtocindex.1 pyshcindextocilm.1 \
	pyshexpanddh.1 pyshexpandglq.1 pyshexpandlsq.1 pyshfindlwin.1 \
	pyshlocalizedadmitcorr.1 pyshmtvaropt.1 pyshmtvaropt0.1 pyshmultiply.1 pywl.1 pywlcurv.1 \
	pyshmultitapercse.1 pyshmultitaperse.1 pyshpowerl.1 pyshpowerdensityl.1 pyshcrosspowerdensityl.1 \
	pyshpowerspectrum.1 pyshpowerspectrumdensity.1 pyshcrosspowerspectrum.1 \
	pyshcrosspowerspectrumdensity.1 pyshcrosspowerl.1 pyshread.1 pyshreaderror.1 pyshreadjpl.1 pyshreadjplerror.1 \
	pyshread2.1 pyshread2error.1 \
	pyshreturntapers.1 pyshreturntapersm.1 pyshrotatecoef.1 pyshrotaterealcoef.1 \
	pywigner3j.1 pyshsjkpg.1 pyshsjkpg0.1 pysphericalcapcoef.1 pyshmtdebias.1 pymakegravgriddh.1 \
	pyshconfidence.1 pyshmagpowerl.1 pyshmagpowerspectrum.1 pymakegriddhc.1 pyshexpanddhc.1 \
	pymakegridglqc.1 pyshexpandglqc.1 pyshpowerlc.1 pyshpowerdensitylc.1 pyshcrosspowerdensitylc.1 \
	pyshpowerspectrumc.1 pyshpowerspectrumdensityc.1 pyshcrosspowerspectrumc.1 pyshcrosspowerspectrumdensityc.1 \
	pyshcrosspowerlc.1 pyshbiasadmitcorr.1 pyshcilmtovector.1 pyshvectortocilm.1 pyyilmindex.1 \
	pycomputedmap.1 pyshreturntapersmap.1 pycurve2mask.1 pymakeellipsecoord.1 pycilmplusrhohdh.1 pycilmplusrhohglq.1 \
	pyhilmrhohdh.1 pyhilmrhohglq.1 pynormalgravity.1 pymakegravgradgriddh.1 pyshreadh.1 pyshreaderrorh.1

HTMLFILES_legendre =  pyplmbar.html
#pyplbar_d1.html pyplbar.html pyplegendre.html pyplegendre_d1.html \
	pyplegendrea_d1.html pyplegendrea.html pyplmbar_d1.html  pyplmindex.html \
	pyplon.html pyplon_d1.html pyplmon.html pyplmon_d1.html pyplschmidt.html pyplschmidt_d1.html \
	pyplmschmidt.html pyplmschmidt_d1.html 

HTMLFILES_transformations = pymakegriddh.html 
#pydhaj.html  pymakegrid2d.html pyglqgridcoord.html \
	 pymakegridglq.html pymakegridpoint.html pyprecompute.html \
	pyshexpanddh.html pyshexpandglq.html pyshexpandlsq.html pyshmultiply.html \
	pyshexpanddhc.html pymakegriddhc.html pymakegridglqc.html pyshexpandglqc.html \
	pypreglq.html pynglq.html pynglqsh.html pynglqshn.html

HTMLFILES_io = 
#pyshread.html pyshreaderror.html pyshreadjpl.html pyshreadjplerror.html pyshread2.html \
 	pyshread2error.html pyshrtoc.html pyshctor.html pyshcilmtocindex.html pyshcindextocilm.html \
 	pyshcilmtovector.html pyshvectortocilm.html pyyilmindex.html pyshreadh.html pyshreaderrorh.html
	
HTMLFILES_gspectra = 
#pyshadmitcorr.html pyshpowerl.html \
	pyshpowerdensityl.html pyshcrosspowerdensityl.html pyshpowerspectrum.html \
	pyshpowerspectrumdensity.html pyshcrosspowerspectrum.html pyshcrosspowerspectrumdensity.html \
	pyshcrosspowerl.html pyshconfidence.html pyshpowerlc.html pyshpowerdensitylc.html \
	pyshcrosspowerdensitylc.html pyshpowerspectrumc.html pyshpowerspectrumdensityc.html \
	pyshcrosspowerspectrumc.html pyshcrosspowerspectrumdensityc.html pyshcrosspowerlc.html

HTMLFILES_localized = 
#pycomputed0.html pycomputedg82.html pycomputedm.html \
	pyshbias.html pyshbiask.html pyshfindlwin.html \
	pyshlocalizedadmitcorr.html pyshmtvaropt.html pyshmtvaropt0.html pyshmultitapercse.html \
	pyshmultitaperse.html pyshreturntapers.html pyshreturntapersm.html pyshsjkpg.html \
	pyshsjkpg0.html pysphericalcapcoef.html pyshmtdebias.html pyshbiasadmitcorr.html \
	pycomputedmap.html pyshreturntapersmap.html pycurve2mask.html

HTMLFILES_rotations = 
#pydjpi2.html pyshrotatecoef.html pyshrotaterealcoef.html

HTMLFILES_gravmag = 
#pycilmplusdh.html pycilmplusglq.html pyhilmdh.html pyhilmglq.html pymakegeoidgriddh.html \
	pymakegeoidgridglq.html pymakegeoidgrid2d.html pymakemaggriddh.html pywl.html pywlcurv.html pymakegravgriddh.html \
	pyshmagpowerl.html pyshmagpowerspectrum.html pycilmplusrhohdh.html pycilmplusrhohglq.html pyhilmrhohdh.html \
	pyhilmrhohglq.html pynormalgravity.html pymakegravgradgriddh.html

HTMLFILES_other =  
#pyeigvalsym.html pyeigvalvecsym.html pyeigvalvecsymtri.html \
	pymakecirclecoord.html pyrandomn.html pyrandomgaussian.html \
	pywigner3j.html pymakeellipsecoord.html
 

all: install-man install-html move


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
	
$(HTMLFILES_legendre) : %.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_legendre.html test.html tail_legendre.html > $@
	@rm -f test.html

$(HTMLFILES_transformations) : %.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_transformations.html test.html tail_transformations.html > $@
	@rm -f test.html
	
$(HTMLFILES_io) : %.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_io.html test.html tail_io.html > $@
	@rm -f test.html

$(HTMLFILES_gspectra) : %.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_gspectra.html test.html tail_gspectra.html > $@
	@rm -f test.html

$(HTMLFILES_localized) : %.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_localized.html test.html tail_localized.html > $@
	@rm -f test.html

$(HTMLFILES_rotations) : %.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_rotations.html test.html tail_rotations.html > $@
	@rm -f test.html

$(HTMLFILES_gravmag) : %.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_gravmag.html test.html tail_gravmag.html > $@
	@rm -f test.html

$(HTMLFILES_other) : %.html : %.md
	$(PANDOC) $(PDHTMLFLAGS) $< -o test.html
	cat header_other.html test.html tail_other.html > $@
	@rm -f test.html


.PHONY: move clean uninstall

move:	
	@mv *.1 $(MANDIR)
	@mv $(HTMLFILES_legendre) $(WWWDIR)
#	@mv $(HTMLFILES_transformations) $(WWWDIR)
#	@mv $(HTMLFILES_io) $(WWWDIR)
#	@mv $(HTMLFILES_gspectra) $(WWWDIR)
#	@mv $(HTMLFILES_localized) $(WWWDIR)
#	@mv $(HTMLFILES_rotations) $(WWWDIR)
#	@mv $(HTMLFILES_gravmag) $(WWWDIR)
#	@mv $(HTMLFILES_other) $(WWWDIR)


clean:
	-rm -f $(MANFILES) $(HTMLFILES) 
	
uninstall:
	-cd /usr/local/share/man/man1/ ; rm $(MANFILES)


.SUFFIXES: $(SUFFIXES)  .md .html .1
	

.md.1:	
	printf "%% $(basename $<)(1) SHTOOLS $(VERSION) | SHTOOLS $(VERSION)\n%%\n\%% %s" > temp.md
	stat -t "%F" -f "%Sm" $< >> temp.md
	cat $< >> temp.md
	$(PANDOC) -s -t man temp.md -o $@ --filter delink.hs
	rm -f temp.md

