######################################################################
#
#	make all
#		Run all Python tests.
#
#	make clean
#		Remove compiled Python files and output.
#
#	Written by Mark Wieczorek (March 2015).
#
#######################################################################

PYTHON = python

SHELL=/bin/sh

.PHONY: all clean

all:
	cd ClassInterface ; $(PYTHON) ./ClassExample.py
	cd ClassInterface ; $(PYTHON) ./WindowExample.py
	cd GlobalSpectralAnalysis ; $(PYTHON) ./SHComplexSpectralAnalysis.py
	cd GlobalSpectralAnalysis ; $(PYTHON) ./SHRealSpectralAnalysis.py
	cd IOStorageConversions ; $(PYTHON) ./SHConversions.py
	cd IOStorageConversions ; $(PYTHON) ./SHStorage.py
	cd LocalizedSpectralAnalysis ; $(PYTHON) ./SHMultitaperSE.py
	cd LocalizedSpectralAnalysis ; $(PYTHON) ./SHWindowsBiasOther.py
	cd SHRotations ; $(PYTHON) ./SHRotations.py
	cd GravMag ; $(PYTHON) ./TestGrav.py
	cd GravMag ; $(PYTHON) ./TestCT.py
	cd Other ; $(PYTHON) ./TestOther.py
	cd TestLegendre ; $(PYTHON) ./TestLegendre.py
	cd TimingAccuracy ; $(PYTHON) ./TimingAccuracyDH.py
	cd TimingAccuracy ; $(PYTHON) ./TimingAccuracyDHC.py
	cd TimingAccuracy ; $(PYTHON) ./TimingAccuracyGLQ.py
	cd TimingAccuracy ; $(PYTHON) ./TimingAccuracyGLQC.py

clean:
	-rm -f ClassInterface/*.pyc
	-rm -f GlobalSpectralAnalysis/*.pyc
	-rm -f GlobalSpectralAnalysis/*.png
	-rm -f IOStorageConversions/*.pyc
	-rm -f IOStorageConversions/*.png
	-rm -f LocalizedSpectralAnalysis/*.pyc
	-rm -f LocalizedSpectralAnalysis/*.png
	-rm -f SHRotations/*.pyc
	-rm -f SHRotations/*.png
	-rm -f TestLegendre/*.pyc
	-rm -f TestLegendre/*.png
	-rm -f TimingAccuracy/*.pyc
	-rm -f TimingAccuracy/*.png
	-rm -f GravMag/*.pyc
	-rm -f GravMag/*.png
	-rm -f Other/*.pyc
	-rm -f Other/*.png
