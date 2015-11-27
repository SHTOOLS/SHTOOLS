################################################################################
#
#	make all
#		Compile example programs and tests in their directories.
#
#	make run-fortran-tests
#		Run the Fortran examples and test suite
#
# 	make clean
#		Delete the compile executibles.
#
#	Written by Mark Wieczorek (June 2012).
#
################################################################################

F95 = gfortran

FFTW = -lfftw3
LAPACK = -llapack
BLAS = -lblas

LIBPATH = ../../lib
MODPATH = ../../modules/

LIBNAME = SHTOOLS

SHELL=/bin/sh
MAKE = make

.PHONY: all getflags clean run-fortran-tests


all: getflags

	@echo
	@echo -------------------------------------------------------------------------------------------------------
	@echo Compiling example and test suite code in SHTOOLS/examples/fortran with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo -------------------------------------------------------------------------------------------------------
	@echo
	
	$(F95) MarsCrustalThickness/MarsCrustalThickness.f95 $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS) -o MarsCrustalThickness/MarsCrustalThickness
	
	$(F95) SHCilmPlus/TestCilmPlus.f95 $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS) -o SHCilmPlus/TestCilmPlus
	
	$(F95) SHExpandDH/TestExpandDH.f95 $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS) -o SHExpandDH/TestExpandDH
	
	$(F95) SHExpandLSQ/TestSHExpandLSQ.f95 $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS) -o SHExpandLSQ/TestSHExpandLSQ
	
	$(F95) SHLocalizedAdmitCorr/LocalizedAdmitCorr.f95 $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS) -o SHLocalizedAdmitCorr/LocalizedAdmitCorr
	
	$(F95) SHMag/SHMag.f95 $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS) -o SHMag/SHMag
	
	$(F95) SHRotate/TestSHRotate.f95 $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS) -o SHRotate/TestSHRotate
	
	$(F95) TimingAccuracy/TimingAccuracyDH.f95 $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS) -o TimingAccuracy/TimingAccuracyDH
	
	$(F95) TimingAccuracy/TimingAccuracyDHC.f95 $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS) -o TimingAccuracy/TimingAccuracyDHC
	
	$(F95) TimingAccuracy/TimingAccuracyGLQ.f95 $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS) -o TimingAccuracy/TimingAccuracyGLQ
	
	$(F95) TimingAccuracy/TimingAccuracyGLQC.f95 $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS) -o TimingAccuracy/TimingAccuracyGLQC
	
run-fortran-tests:
	 @echo MARS-CRUSTAL-THICKNESS
	 cd MarsCrustalThickness ; ./MarsCrustalThickness < input.txt
	 @echo
	 @echo TEST-CILMPLUS
	 cd SHCilmPlus ; ./TestCilmPlus
	 @echo
	 @echo TEST-EXPANDDH
	 cd SHExpandDH ; ./TestExpandDH
	 @echo
	 @echo TEST-SHEXPANDLSQ
	 cd SHExpandLSQ ; ./TestSHExpandLSQ < input.txt
	 @echo
	 @echo LOCALIZED-ADMIT-CORR
	 cd SHLocalizedAdmitCorr ; ./LocalizedAdmitCorr < input.txt
	 @echo
	 @echo SHMAG
	 cd SHMag ; ./SHMag
	 @echo
	 @echo TEST-SHROTATE
	 cd SHRotate ; ./TestSHRotate < input.txt
	 @echo
	 @echo "*** The following 8 timing and accuracy tests will take about an hour to complete ***"
	 @echo "*** If you are impatient, press Control-C to exit ***"
	 @echo
	 @echo TIMING-ACCURACY-GLQ
	 cd TimingAccuracy ; ./TimingAccuracyGLQ < input1 > output1
	 @echo
	 @echo TIMING-ACCURACY-GLQ
	 cd TimingAccuracy ; ./TimingAccuracyGLQ < input2 > output2
	 @echo
	 @echo TIMING-ACCURACY-DH
	 cd TimingAccuracy ; ./TimingAccuracyDH < input3 > output3
	 @echo
	 @echo TIMING-ACCURACY-DH
	 cd TimingAccuracy ; ./TimingAccuracyDH < input4 > output4
	 @echo
	 @echo TIMING-ACCURACY-DHC
	 cd TimingAccuracy ; ./TimingAccuracyDHC < input5 > output5
	 @echo
	 @echo TIMING-ACCURACY-DHC
	 cd TimingAccuracy ; ./TimingAccuracyDHC < input6 > output6
	 @echo
	 @echo TIMING-ACCURACY-GLQC
	 cd TimingAccuracy ; ./TimingAccuracyGLQC < input7 > output7
	 @echo
	 @echo TIMING-ACCURACY-GLQC
	 cd TimingAccuracy ; ./TimingAccuracyGLQC < input8 > output8

getflags:

ifeq ($(F95),f95)
F95FLAGS ?= -m64 -O3 -YEXT_NAMES=LCS -YEXT_SFX=_ -fpic
MODFLAG = -p $(MODPATH)
endif

ifeq ($(F95),gfortran)
F95FLAGS ?= -m64 -fPIC -O3
MODFLAG = -I$(MODPATH)
endif

ifeq ($(F95),ifort)
F95FLAGS ?= -m64 -free -O3 -Tf
MODFLAG = -I$(MODPATH)
endif

ifeq ($(F95),g95)
F95FLAGS ?= -O3 -fno-second-underscore 
MODFLAG = -I$(MODPATH)
endif

ifeq ($(F95),pgf90)
F95FLAGS ?= -fast 
MODFLAG = -Imodpath
endif

ifeq ($(origin F95FLAGS), undefined)
F95FLAGS = -m64 -O3
MODFLAG = -I$(MODPATH)
endif


clean:
	-rm -f MarsCrustalThickness/MarsCrustalThickness
	-rm -f MarsCrustalThickness/*.sh
	-rm -f MarsCrustalThickness/*.dat
	-rm -f SHCilmPlus/TestCilmPlus
	-rm -f SHCilmPlus/*.dat
	-rm -f SHExpandDH/TestExpandDH
	-rm -f SHExpandLSQ/TestSHExpandLSQ
	-rm -f SHLocalizedAdmitCorr/LocalizedAdmitCorr
	-rm -f SHLocalizedAdmitCorr/*.dat
	-rm -f SHMag/SHMag
	-rm -f SHMag/*.dat
	-rm -f SHRotate/TestSHRotate
	-rm -f SHRotate/*.sh
	-rm -f SHRotate/*.dat
	-rm -f TimingAccuracy/TimingAccuracyDH
	-rm -f TimingAccuracy/TimingAccuracyDHC
	-rm -f TimingAccuracy/TimingAccuracyGLQ
	-rm -f TimingAccuracy/TimingAccuracyGLQC
	-rm -f TimingAccuracy/output*
	-rm -f TimingAccuracy/*.maxerror
	-rm -f TimingAccuracy/*.rmserror
	-rm -f TimingAccuracy/*.timef
	-rm -f TimingAccuracy/*.timei
	


