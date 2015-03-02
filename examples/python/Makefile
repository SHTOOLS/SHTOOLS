######################################################################
#
#	make all
#		Run all Python tests.
#
#	make clean
#		Remove compiled Python files.
#
#	Written by Mark Wieczorek (March 2015).
#
#######################################################################

PYTHON = python

SHELL=/bin/tcsh

MAKE = make

.PHONY: all clean


all:
	$(PYTHON) GlobalSpectralAnalysis/SHComplexSpectralAnalysis.py
	$(PYTHON) GlobalSpectralAnalysis/SHRealSpectralAnalysis.py
	$(PYTHON) IOStorageConversions/SHConversions.py
	$(PYTHON) IOStorageConversions/SHStorage.py
	$(PYTHON) LocalizedSpectralAnalysis/SHMultitaperSE.py
	$(PYTHON) LocalizedSpectralAnalysis/SHWindowsBiasOther.py
	$(PYTHON) SHRotations/SHRotations.py
	$(PYTHON) TestLegendre/TestLegendre.py
	$(PYTHON) TimingAccuracy/TimingAccuracyDH.py
	$(PYTHON) TimingAccuracy/TimingAccuracyDHC.py
	$(PYTHON) TimingAccuracy/TimingAccuracyGLQ.py
	$(PYTHON) TimingAccuracy/TimingAccuracyGLQC.py

clean:
	-rm -f GlobalSpectralAnalysis/*.pyc
	-rm -f IOStorageConversions/*.pyc
	-rm -f LocalizedSpectralAnalysis/*.pyc
	-rm -f SHRotations/*.pyc
	-rm -f TestLegendre/*.pyc
	-rm -f TimingAccuracy/*.pyc


