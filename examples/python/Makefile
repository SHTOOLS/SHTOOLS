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
	cd GlobalSpectralAnalysis ; $(PYTHON) ./SHComplexSpectralAnalysis.py
	cd GlobalSpectralAnalysis ; $(PYTHON) ./SHRealSpectralAnalysis.py
	cd IOStorageConversions ; $(PYTHON) ./SHConversions.py
	cd IOStorageConversions ; $(PYTHON) ./SHStorage.py
	cd LocalizedSpectralAnalysis ; $(PYTHON) ./SHMultitaperSE.py
	cd LocalizedSpectralAnalysis ; $(PYTHON) ./SHWindowsBiasOther.py
	cd SHRotations ; $(PYTHON) ./SHRotations.py
	cd TestLegendre ; $(PYTHON) ./TestLegendre.py
	cd TimingAccuracy ; $(PYTHON) ./TimingAccuracyDH.py
	cd TimingAccuracy ; $(PYTHON) ./TimingAccuracyDHC.py
	cd TimingAccuracy ; $(PYTHON) ./TimingAccuracyGLQ.py
	cd TimingAccuracy ; $(PYTHON) ./TimingAccuracyGLQC.py

clean:
	-rm -f GlobalSpectralAnalysis/*.pyc
	-rm -f IOStorageConversions/*.pyc
	-rm -f LocalizedSpectralAnalysis/*.pyc
	-rm -f SHRotations/*.pyc
	-rm -f TestLegendre/*.pyc
	-rm -f TimingAccuracy/*.pyc


