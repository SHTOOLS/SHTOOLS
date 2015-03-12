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

MAN2HTML = man2html
MAN2HTMLFLAGS = -compress -headmap headmapfile -uelem B -topm 0 -bare -nodepage 
#-cgiurl "file:$${title}.$${section}$${subsection}.html"
POD2MAN = pod2man
POD2MANFLAGS = -c "SHTOOLS $(VERSION)" -r "SHTOOLS $(VERSION)"
MANDIR = ../../man/man1
WWWDIR = ../../www/man/python
MANWIDTH = 130

PODFILES = cilmplusdh.pod cilmplusglq.pod computed0.pod computedg82.pod computedm.pod \
	dhaj.pod djpi2.pod eigvalsym.pod eigvalvecsym.pod eigvalvecsymtri.pod \
	hilmdh.pod hilmglq.pod makecirclecoord.pod makegeoidgriddh.pod makegeoidgridglq.pod \
	makegeoidgrid2d.pod makegrid2d.pod glqgridcoord.pod \
	makegriddh.pod makegridglq.pod makegridpoint.pod makemaggriddh.pod \
	planetsconstants.pod plbar_d1.pod plbar.pod plegendre.pod plegendre_d1.pod \
	plegendrea_d1.pod plegendrea.pod plmbar_d1.pod plmbar.pod plmindex.pod \
	plon.pod plon_d1.pod plmon.pod plmon_d1.pod plschmidt.pod plschmidt_d1.pod \
	plmschmidt.pod plmschmidt_d1.pod precompute.pod preglq.pod nglq.pod nglqsh.pod \
	nglqshn.pod randomn.pod randomgaussian.pod \
	shadmitcorr.pod shbias.pod shbiask.pod shrtoc.pod \
	shctor.pod shcilmtocindex.pod shcindextocilm.pod  \
	shexpanddh.pod shexpandglq.pod shexpandlsq.pod shfindlwin.pod \
	shlocalizedadmitcorr.pod shmtvaropt.pod shmtvaropt0.pod shmultiply.pod \
	wl.pod wl_curv.pod shmultitapercse.pod shmultitaperse.pod shpowerl.pod \
	shpowerdensityl.pod shcrosspowerdensityl.pod shpowerspectrum.pod \
	shpowerspectrumdensity.pod shcrosspowerspectrum.pod shcrosspowerspectrumdensity.pod \
	shcrosspowerl.pod shread.pod shreaderror.pod shreadjpl.pod shreadjplerror.pod shread2.pod \
	shread2error.pod shreturntapers.pod shreturntapersm.pod shrotatecoef.pod shrotaterealcoef.pod \
	wigner3j.pod shsjkpg.pod shsjkpg0.pod shtools.pod sphericalcapcoef.pod shmtdebias.pod \
	makegravgriddh.pod shconfidence.pod shmagpowerl.pod shmagpowerspectrum.pod \
	makegriddhc.pod shexpanddhc.pod makegridglqc.pod shexpandglqc.pod shpowerlc.pod \
	shpowerdensitylc.pod shcrosspowerdensitylc.pod shpowerspectrumc.pod \
	shpowerspectrumdensityc.pod shcrosspowerspectrumc.pod shcrosspowerspectrumdensityc.pod \
	shcrosspowerlc.pod shbiasadmitcorr.pod shcilmtovector.pod shvectortocilm.pod \
	yilmindex.pod computedmap.pod shreturntapers.pod curve2mask.pod makeellipsecoord.pod \
	cilmplusrhohdh.pod cilmplusglq.pod hilmrhohdh.pod hilmrhohglq.pod normalgravity.pod \
	makegravgradgriddh.pod shreadh.pod shreaderrorh.pod

MANFILES = cilmplusdh.1 cilmplusglq.1 computed0.1 computedg82.1 computedm.1 \
	dhaj.1 djpi2.1 eigvalsym.1 eigvalvecsym.1 eigvalvecsymtri.1 hilmdh.1 hilmglq.1 \
	makecirclecoord.1 makegeoidgriddh.1 makegeoidgridglq.1 makegeoidgrid2d.1 \
	makegrid2d.1 glqgridcoord.1 makegriddh.1 \
	makegridglq.1 makegridpoint.1 makemaggriddh.1 planetsconstants.1 plbar_d1.1 \
	plbar.1 plegendre.1 plegendre_d1.1 plegendrea_d1.1 plegendrea.1 plmbar_d1.1 \
	plmbar.1 plmindex.1 plon.1 plon_d1.1 plmon.1 plmon_d1.1 plschmidt.1 plschmidt_d1.1 \
	plmschmidt.1 plmschmidt_d1.1 precompute.1 preglq.1 nglq.1 nglqsh.1 nglqshn.1 \
	randomn.1 randomgaussian.1 \
	shadmitcorr.1 shbias.1 shbiask.1 shrtoc.1 shctor.1 shcilmtocindex.1 shcindextocilm.1 \
	shexpanddh.1 shexpandglq.1 shexpandlsq.1 shfindlwin.1 \
	shlocalizedadmitcorr.1 shmtvaropt.1 shmtvaropt0.1 shmultiply.1 wl.1 wlcurv.1 \
	shmultitapercse.1 shmultitaperse.1 shpowerl.1 shpowerdensityl.1 shcrosspowerdensityl.1 \
	shpowerspectrum.1 shpowerspectrumdensity.1 shcrosspowerspectrum.1 \
	shcrosspowerspectrumdensity.1 shcrosspowerl.1 shread.1 shreaderror.1 shreadjpl.1 shreadjplerror.1 \
	shread2.1 shread2error.1 \
	shreturntapers.1 shreturntapersm.1 shrotatecoef.1 shrotaterealcoef.1 \
	wigner3j.1 shsjkpg.1 shsjkpg0.1 shtools.1 sphericalcapcoef.1 shmtdebias.1 makegravgriddh.1 \
	shconfidence.1 shmagpowerl.1 shmagpowerspectrum.1 makegriddhc.1 shexpanddhc.1 \
	makegridglqc.1 shexpandglqc.1 shpowerlc.1 shpowerdensitylc.1 shcrosspowerdensitylc.1 \
	shpowerspectrumc.1 shpowerspectrumdensityc.1 shcrosspowerspectrumc.1 shcrosspowerspectrumdensityc.1 \
	shcrosspowerlc.1 shbiasadmitcorr.1 shcilmtovector.1 shvectortocilm.1 yilmindex.1 \
	computedmap.1 shreturntapersmap.1 curve2mask.1 makeellipsecoord.1 cilmplusrhohdh.1 cilmplusrhohglq.1 \
	hilmrhohdh.1 hilmrhohglq.1 normalgravity.1 makegravgradgriddh.1 shreadh.1 shreaderrorh.1

PYMANFILES = pycilmplusdh.1 pycilmplusglq.1 pycomputed0.1 pycomputedg82.1 pycomputedm.1 \
	pydhaj.1 pydjpi2.1 pyeigvalsym.1 pyeigvalvecsym.1 pyeigvalvecsymtri.1 pyhilmdh.1 pyhilmglq.1 \
	pymakecirclecoord.1 pymakegeoidgriddh.1 pymakegeoidgridglq.1 pymakegeoidgrid2d.1 \
	pymakegrid2d.1 pyglqgridcoord.1 pymakegriddh.1 \
	pymakegridglq.1 pymakegridpoint.1 pymakemaggriddh.1 pyplanetsconstants.1 pyplbar_d1.1 \
	pyplbar.1 pyplegendre.1 pyplegendre_d1.1 pyplegendrea_d1.1 pyplegendrea.1 pyplmbar_d1.1 \
	pyplmbar.1 pyplmindex.1 pyplon.1 pyplon_d1.1 pyplmon.1 pyplmon_d1.1 pyplschmidt.1 pyplschmidt_d1.1 \
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
	pywigner3j.1 pyshsjkpg.1 pyshsjkpg0.1 pyshtools.1 pysphericalcapcoef.1 pyshmtdebias.1 pymakegravgriddh.1 \
	pyshconfidence.1 pyshmagpowerl.1 pyshmagpowerspectrum.1 pymakegriddhc.1 pyshexpanddhc.1 \
	pymakegridglqc.1 pyshexpandglqc.1 pyshpowerlc.1 pyshpowerdensitylc.1 pyshcrosspowerdensitylc.1 \
	pyshpowerspectrumc.1 pyshpowerspectrumdensityc.1 pyshcrosspowerspectrumc.1 pyshcrosspowerspectrumdensityc.1 \
	pyshcrosspowerlc.1 pyshbiasadmitcorr.1 pyshcilmtovector.1 pyshvectortocilm.1 pyyilmindex.1 \
	pycomputedmap.1 pyshreturntapersmap.1 pycurve2mask.1 pymakeellipsecoord.1 pycilmplusrhohdh.1 pycilmplusrhohglq.1 \
	pyhilmrhohdh.1 pyhilmrhohglq.1 pynormalgravity.1 pymakegravgradgriddh.1 pyshreadh.1 pyshreaderrorh.1

HTMLFILES_modules = planetsconstants.html shtools.html

HTMLFILES_legendre =  plbar_d1.html plbar.html plegendre.html plegendre_d1.html \
	plegendrea_d1.html plegendrea.html plmbar_d1.html plmbar.html plmindex.html \
	plon.html plon_d1.html plmon.html plmon_d1.html plschmidt.html plschmidt_d1.html \
	plmschmidt.html plmschmidt_d1.html 

HTMLFILES_transformations = dhaj.html  makegrid2d.html glqgridcoord.html \
	makegriddh.html makegridglq.html makegridpoint.html precompute.html \
	shexpanddh.html shexpandglq.html shexpandlsq.html shmultiply.html \
	shexpanddhc.html makegriddhc.html makegridglqc.html shexpandglqc.html \
	preglq.html nglq.html nglqsh.html nglqshn.html

HTMLFILES_io = shread.html shreaderror.html shreadjpl.html shreadjplerror.html shread2.html \
 	shread2error.html shrtoc.html shctor.html shcilmtocindex.html shcindextocilm.html \
 	shcilmtovector.html shvectortocilm.html yilmindex.html shreadh.html shreaderrorh.html
	
HTMLFILES_gspectra = shadmitcorr.html shpowerl.html \
	shpowerdensityl.html shcrosspowerdensityl.html shpowerspectrum.html \
	shpowerspectrumdensity.html shcrosspowerspectrum.html shcrosspowerspectrumdensity.html \
	shcrosspowerl.html shconfidence.html shpowerlc.html shpowerdensitylc.html \
	shcrosspowerdensitylc.html shpowerspectrumc.html shpowerspectrumdensityc.html \
	shcrosspowerspectrumc.html shcrosspowerspectrumdensityc.html shcrosspowerlc.html

HTMLFILES_localized = computed0.html computedg82.html computedm.html \
	shbias.html shbiask.html shfindlwin.html \
	shlocalizedadmitcorr.html shmtvaropt.html shmtvaropt0.html shmultitapercse.html \
	shmultitaperse.html shreturntapers.html shreturntapersm.html shsjkpg.html \
	shsjkpg0.html sphericalcapcoef.html shmtdebias.html shbiasadmitcorr.html \
	computedmap.html shreturntapersmap.html curve2mask.html

HTMLFILES_rotations = djpi2.html shrotatecoef.html shrotaterealcoef.html

HTMLFILES_gravmag = cilmplusdh.html cilmplusglq.html hilmdh.html hilmglq.html makegeoidgriddh.html \
	makegeoidgridglq.html makegeoidgrid2d.html makemaggriddh.html wl.html wlcurv.html makegravgriddh.html \
	shmagpowerl.html shmagpowerspectrum.html cilmplusrhohdh.html cilmplusrhohglq.html hilmrhohdh.html \
	hilmrhohglq.html normalgravity.html makegravgradgriddh.html

HTMLFILES_other =  eigvalsym.html eigvalvecsym.html eigvalvecsymtri.html \
	makecirclecoord.html randomn.html randomgaussian.html \
	wigner3j.html makeellipsecoord.html
 

all: install-man install-html move


install-man: $(MANFILES)
	@echo
	@echo MAN PAGES SUCCESSFULLY MADE
	@echo
	

install-html: $(HTMLFILES_modules) $(HTMLFILES_legendre) $(HTMLFILES_transformations) \
	$(HTMLFILES_io) $(HTMLFILES_gspectra) $(HTMLFILES_localized) \
	$(HTMLFILES_rotations) $(HTMLFILES_gravmag) \
	$(HTMLFILES_other)
	@echo
	@echo HTML WEB PAGES SUCCESSFULLY MADE
	@echo
	

$(HTMLFILES_modules) : %.html : py%.1
	@rm -f test.html
	man ./$< | $(MAN2HTML) $(MAN2HTMLFLAGS) > test.html
	cat header_modules.html test.html tail_modules.html > $@
	@rm -f test.html
	
$(HTMLFILES_legendre) : %.html : py%.1
	@rm -f test.html
	man ./$< | $(MAN2HTML) $(MAN2HTMLFLAGS) > test.html
	cat header_legendre.html test.html tail_legendre.html > $@
	@rm -f test.html

$(HTMLFILES_transformations) : %.html : py%.1
	@rm -f test.html
	man ./$< | $(MAN2HTML) $(MAN2HTMLFLAGS) > test.html
	cat header_transformations.html test.html tail_transformations.html > $@
	@rm -f test.html
	
$(HTMLFILES_io) : %.html : py%.1
	@rm -f test.html
	man ./$< | $(MAN2HTML) $(MAN2HTMLFLAGS) > test.html
	cat header_io.html test.html tail_io.html > $@
	@rm -f test.html

$(HTMLFILES_gspectra) : %.html : py%.1
	@rm -f test.html
	man ./$< | $(MAN2HTML) $(MAN2HTMLFLAGS) > test.html
	cat header_gspectra.html test.html tail_gspectra.html > $@
	@rm -f test.html

$(HTMLFILES_localized) : %.html : py%.1
	@rm -f test.html
	man ./$< | $(MAN2HTML) $(MAN2HTMLFLAGS) > test.html
	cat header_localized.html test.html tail_localized.html > $@
	@rm -f test.html

$(HTMLFILES_rotations) : %.html : py%.1
	@rm -f test.html
	man ./$< | $(MAN2HTML) $(MAN2HTMLFLAGS) > test.html
	cat header_rotations.html test.html tail_rotations.html > $@
	@rm -f test.html

$(HTMLFILES_gravmag) : %.html : py%.1
	@rm -f test.html
	man ./$< | $(MAN2HTML) $(MAN2HTMLFLAGS) > test.html
	cat header_gravmag.html test.html tail_gravmag.html > $@
	@rm -f test.html

$(HTMLFILES_other) : %.html : py%.1
	@rm -f test.html
	man ./$< | $(MAN2HTML) $(MAN2HTMLFLAGS) > test.html
	cat header_other.html test.html tail_other.html > $@
	@rm -f test.html


.PHONY: move clean uninstall

move:	
	@mv *.1 $(MANDIR)
	@mv $(HTMLFILES_modules) $(WWWDIR)
	@mv $(HTMLFILES_legendre) $(WWWDIR)
	@mv $(HTMLFILES_transformations) $(WWWDIR)
	@mv $(HTMLFILES_io) $(WWWDIR)
	@mv $(HTMLFILES_gspectra) $(WWWDIR)
	@mv $(HTMLFILES_localized) $(WWWDIR)
	@mv $(HTMLFILES_rotations) $(WWWDIR)
	@mv $(HTMLFILES_gravmag) $(WWWDIR)
	@mv $(HTMLFILES_other) $(WWWDIR)


clean:
	-rm -f $(MANFILES) $(HTMLFILES) 
	
uninstall:
	-cd /usr/local/share/man/man1/ ; rm $(PYMANFILES)


.SUFFIXES: $(SUFFIXES)  .pod .html .1
	

.pod.1:
	$(POD2MAN) $(POD2MANFLAGS) $< > py$@

