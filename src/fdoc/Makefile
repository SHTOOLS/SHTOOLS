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

MDFILES = shtools.md planetsconstants.md plmbar.md plmbar_d1.md plbar.md \
	plbar_d1.md plon.md plmon.md plmon_d1.md plon_d1.md plmschmidt.md \
	plschmidt.md plmschmidt_d1.md plschmidt_d1.md plegendrea.md \
	plegendrea_d1.md plegendre.md plegendre_d1.md plmindex.md shexpanddh.md \
	makegriddh.md shexpanddhc.md makegriddhc.md shglq.md glqgridcoord.md \
	shexpandglq.md makegridglq.md shexpandglqc.md makegridglqc.md \
	shexpandlsq.md makegrid2d.md makegridpoint.md makegridpointc.md \
	shmultiply.md shread.md  shread2.md shreadjpl.md shcilmtocindex.md \
	shcindextocilm.md shcilmtovector.md shvectortocilm.md yilmindexvector.md \
	shrtoc.md shctor.md shpowerl.md shpowerdensityl.md shcrosspowerl.md \
	shcrosspowerdensityl.md shpowerspectrum.md shpowerspectrumdensity.md \
	shcrosspowerspectrum.md shcrosspowerspectrumdensity.md shadmitcorr.md \
	shconfidence.md shpowerlc.md shpowerdensitylc.md shcrosspowerlc.md \
	shcrosspowerdensitylc.md shpowerspectrumc.md shpowerspectrumdensityc.md \
	shcrosspowerspectrumc.md shcrosspowerspectrumdensityc.md shmultitaperse.md \
	shmultitapercse.md shlocalizedadmitcorr.md shreturntapers.md \
	shreturntapersm.md computedm.md computedg82.md shfindlwin.md shbiask.md \
	shbiasadmitcorr.md shmtdebias.md shmtvaropt.md shsjkpg.md \
	shreturntapersmap.md computedmap.md curve2mask.md shbias.md \
	sphericalcapcoef.md djpi2.md shrotatecoef.md shrotaterealcoef.md \
	makecirclecoord.md makeellipsecoord.md randomn.pod randomgaussian.md \
	preglq.md eigvalvecsym.md eigvalvecsymtri.md eigvalsym.md wigner3j.md \
	makegravgriddh.md makegravgradgriddh.md makegeoidgrid.md cilmplus.md \
	cilmminus.md cilmplusrhoh.md cilmminusrhoh.md batohilm.md batohilmrhoh.md \
	downcontfilterma.md downcontfiltermc.md normalgravity.md makemaggriddh.md \
	shmagpowerspectrum.md shmagpowerl.md shmtcouplingmatrix.md \
	shmultitapermaskcse.md shmultitapermaskse.md shbiaskmask.md dhaj.md

MANFILES = $(addprefix $(MANDIR)/, shtools.1 planetsconstants.1 plmbar.1 \
	plmbar_d1.1 plbar.1 plbar_d1.1 plon.1 plmon.1 plmon_d1.1 plon_d1.1 \
	plmschmidt.1 plschmidt.1 plmschmidt_d1.1 plschmidt_d1.1 plegendrea.1 \
	plegendrea_d1.1 plegendre.1 plegendre_d1.1 plmindex.1 shexpanddh.1 \
	makegriddh.1 shexpanddhc.1 makegriddhc.1 glqgridcoord.1 shexpandglq.1 \
	makegridglq.1 shexpandglqc.1 makegridglqc.1 shexpandlsq.1 makegrid2d.1 \
	makegridpoint.1 makegridpointc.1 shmultiply.1 shread.1 shread2.1 \
	shreadjpl.1 shcilmtocindex.1 shcindextocilm.1 shcilmtovector.1 \
	shvectortocilm.1 yilmindexvector.1 shrtoc.1 shctor.1 shpowerl.1 \
	shpowerdensityl.1 shcrosspowerl.1 shcrosspowerdensityl.1 \
	shpowerspectrum.1 shpowerspectrumdensity.1 shcrosspowerspectrum.1 \
	shcrosspowerspectrumdensity.1 shadmitcorr.1 shconfidence.1 shpowerlc.1 \
	shpowerdensitylc.1 shcrosspowerlc.1 shcrosspowerdensitylc.1 \
	shpowerspectrumc.1 shpowerspectrumdensityc.1 shcrosspowerspectrumc.1 \
	shcrosspowerspectrumdensityc.1 shmultitaperse.1 shmultitapercse.1 \
	shlocalizedadmitcorr.1 shreturntapers.1 shreturntapersm.1 computedm.1 \
	computedg82.1 shfindlwin.1 shbiask.1 shbiasadmitcorr.1 shmtdebias.1 \
	shmtvaropt.1 shsjkpg.1 shreturntapersmap.1 computedmap.1 curve2mask.1 \
	shbias.1 sphericalcapcoef.1 djpi2.1 shrotatecoef.1 shrotaterealcoef.1 \
	makecirclecoord.1 makeellipsecoord.1 randomn.1 randomgaussian.1 preglq.1 \
	eigvalvecsym.1 eigvalvecsymtri.1 eigvalsym.1 wigner3j.1 makegravgriddh.1 \
	makegravgradgriddh.1 makegeoidgrid.1 cilmplus.1 cilmminus.1 \
	cilmplusrhoh.1 cilmminusrhoh.1 batohilm.1 batohilmrhoh.1 \
	downcontfilterma.1 downcontfiltermc.1 normalgravity.1 makemaggriddh.1 \
	shmagpowerspectrum.1 shmagpowerl.1 shmtcouplingmatrix.1 \
	shmultitapermaskcse.1 shmultitapermaskse.1 shbiaskmask.1 dhaj.1 shglq.1)


all: install-man
	$(PYTHON) make_www_doc.py

install-man: $(MANFILES)
	@echo
	@echo "--> Fortran 95 man pages created successfully"
	@echo


.PHONY: clean uninstall

clean:
	-rm -f $(MANFILES)

uninstall:
	-cd /usr/local/share/man/man1/ ; rm $(MANFILES)


.SUFFIXES: $(SUFFIXES) .md .1

$(MANDIR)/%.1 : %.md
	printf "%% $(basename $<)(1) Fortran 95 | SHTOOLS $(VERSION)\n%%\n%% %s" > temp.md
	stat -t "%F" -f "%Sm" $< >> temp.md
	cat $< >> temp.md
	$(PANDOC) -s -t man temp.md -o $@ --filter delink.hs
	rm -f temp.md
