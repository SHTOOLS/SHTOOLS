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
PREFIX=/usr/local

PANDOC = pandoc
MANDIR = ../../man/man3
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
	shmultitapermaskcse.md shmultitapermaskse.md shbiaskmask.md dhaj.md \
	makemaggradgriddh.md shrotatetapers.md slepiancoeffs.md \
	slepiancoeffstosh.md shscouplingmatrix.md shslepianvar.md \
	shscouplingmatrixcap.md makegravgridpoint.md makemaggridpoint.md \
	makegradientdh.md lsq_g.md

MANFILES = $(addprefix $(MANDIR)/, shtools.3 planetsconstants.3 plmbar.3 \
	plmbar_d1.3 plbar.3 plbar_d1.3 plon.3 plmon.3 plmon_d1.3 plon_d1.3 \
	plmschmidt.3 plschmidt.3 plmschmidt_d1.3 plschmidt_d1.3 plegendrea.3 \
	plegendrea_d1.3 plegendre.3 plegendre_d1.3 plmindex.3 shexpanddh.3 \
	makegriddh.3 shexpanddhc.3 makegriddhc.3 glqgridcoord.3 shexpandglq.3 \
	makegridglq.3 shexpandglqc.3 makegridglqc.3 shexpandlsq.3 makegrid2d.3 \
	makegridpoint.3 makegridpointc.3 shmultiply.3 shread.3 shread2.3 \
	shreadjpl.3 shcilmtocindex.3 shcindextocilm.3 shcilmtovector.3 \
	shvectortocilm.3 yilmindexvector.3 shrtoc.3 shctor.3 shpowerl.3 \
	shpowerdensityl.3 shcrosspowerl.3 shcrosspowerdensityl.3 \
	shpowerspectrum.3 shpowerspectrumdensity.3 shcrosspowerspectrum.3 \
	shcrosspowerspectrumdensity.3 shadmitcorr.3 shconfidence.3 shpowerlc.3 \
	shpowerdensitylc.3 shcrosspowerlc.3 shcrosspowerdensitylc.3 \
	shpowerspectrumc.3 shpowerspectrumdensityc.3 shcrosspowerspectrumc.3 \
	shcrosspowerspectrumdensityc.3 shmultitaperse.3 shmultitapercse.3 \
	shlocalizedadmitcorr.3 shreturntapers.3 shreturntapersm.3 computedm.3 \
	computedg82.3 shfindlwin.3 shbiask.3 shbiasadmitcorr.3 shmtdebias.3 \
	shmtvaropt.3 shsjkpg.3 shreturntapersmap.3 computedmap.3 curve2mask.3 \
	shbias.3 sphericalcapcoef.3 djpi2.3 shrotatecoef.3 shrotaterealcoef.3 \
	makecirclecoord.3 makeellipsecoord.3 randomn.3 randomgaussian.3 preglq.3 \
	eigvalvecsym.3 eigvalvecsymtri.3 eigvalsym.3 wigner3j.3 makegravgriddh.3 \
	makegravgradgriddh.3 makegeoidgrid.3 cilmplus.3 cilmminus.3 \
	cilmplusrhoh.3 cilmminusrhoh.3 batohilm.3 batohilmrhoh.3 \
	downcontfilterma.3 downcontfiltermc.3 normalgravity.3 makemaggriddh.3 \
	shmagpowerspectrum.3 shmagpowerl.3 shmtcouplingmatrix.3 \
	shmultitapermaskcse.3 shmultitapermaskse.3 shbiaskmask.3 dhaj.3 shglq.3 \
	makemaggradgriddh.3 shrotatetapers.3 slepiancoeffs.3 slepiancoeffstosh.3 \
	shscouplingmatrix.3 shslepianvar.3 shscouplingmatrixcap.3 \
	makegravgridpoint.3 makemaggridpoint.3 makegradientdh.3 lsq_g.3)


all: install-man
	$(PYTHON) make_www_doc.py
	@echo "Fortran 95 man pages created successfully"

install-man: $(MANFILES)

.PHONY: clean uninstall

clean:
	-rm -f $(MANFILES)

uninstall:
	-cd $(PREFIX)/share/man/man3/ ; rm $(MANFILES)


.SUFFIXES: $(SUFFIXES) .md .3

$(MANDIR)/%.3 : %.md
	printf "%% $(basename $<)(1) Fortran 95 | SHTOOLS $(VERSION)\n%%\n%% %s" > temp.md
	stat -t "%F" -f "%Sm" $< >> temp.md
	cat $< >> temp.md
	$(PANDOC) -s -t man temp.md -o $@ --filter delink.hs
	rm -f temp.md
