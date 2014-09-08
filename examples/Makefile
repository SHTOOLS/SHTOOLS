######################################################################
#
#	make all
#		Compile example programs in their directories.
#
# 	make remove-examples
#		Delete the compile executibles.
#
#	make run-examples
#
#	Written by Mark Wieczorek (June 2012).
#
#######################################################################

F95 = f95

SHELL=/bin/tcsh

MAKE = make

.PHONY: all getflags clean


all: getflags

	@echo
	@echo -------------------------------------------------------------------------------------------------------
	@echo Compiling example code in SHTOOLS/examples with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm $(F95FLAGS)
	@echo -------------------------------------------------------------------------------------------------------
	@echo
	
	$(F95) fftw_wisdom/create_wisdom.f95 $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm $(F95FLAGS) -o fftw_wisdom/create_wisdom
	
	$(F95) MarsCrustalThickness/MarsCrustalThickness.f95 $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm $(F95FLAGS) -o MarsCrustalThickness/MarsCrustalThickness
	
	$(F95) SHCilmPlus/TestCilmPlus.f95 $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm $(F95FLAGS) -o SHCilmPlus/TestCilmPlus
	
	$(F95) SHExpandDH/TestExpandDH.f95 $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm $(F95FLAGS) -o SHExpandDH/TestExpandDH
	
	$(F95) SHExpandLSQ/TestSHExpandLSQ.f95 $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -llapack -lblas -lfftw3 -lm $(F95FLAGS) -o SHExpandLSQ/TestSHExpandLSQ
	
	$(F95) SHLocalizedAdmitCorr/LocalizedAdmitCorr.f95 $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -llapack -lblas -lfftw3 -lm $(F95FLAGS) -o SHLocalizedAdmitCorr/LocalizedAdmitCorr
	
	$(F95) SHMag/SHMag.f95 $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm $(F95FLAGS) -o SHMag/SHMag
	
	$(F95) SHRotate/TestSHRotate.f95 $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm $(F95FLAGS) -o SHRotate/TestSHRotate
	
	$(F95) TimingAccuracy/TimingAccuracyDH.f95 $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm $(F95FLAGS) -o TimingAccuracy/TimingAccuracyDH
	
	$(F95) TimingAccuracy/TimingAccuracyDHC.f95 $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm $(F95FLAGS) -o TimingAccuracy/TimingAccuracyDHC
	
	$(F95) TimingAccuracy/TimingAccuracyGLQ.f95 $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm $(F95FLAGS) -o TimingAccuracy/TimingAccuracyGLQ
	
	$(F95) TimingAccuracy/TimingAccuracyGLQC.f95 $(MODFLAG) -L../lib -lSHTOOLS$(VERSION) -lfftw3 -lm $(F95FLAGS) -o TimingAccuracy/TimingAccuracyGLQC
	
make run-examples:
	 @echo CREATE_WISDOM
	 cd fftw_wisdom ; ./create_wisdom < input.txt
	 @echo
	 @echo MARS CRUSTAL THICKNESS
	 cd MarsCrustalThickness ; ./MarsCrustalThickness < input.txt
	 @echo
	 @echo TEST CILMPLUS
	 cd SHCilmPlus ; ./TestCilmPlus
	 @echo
	 @echo TEST EXPANDDH
	 cd SHExpandDH ; ./TestExpandDH
	 @echo
	 @echo TEST SHEXPANDLSQ
	 cd SHExpandLSQ ; ./TestSHExpandLSQ < input.txt
	 @echo
	 @echo LOCALIZEDADMITCORR
	 cd SHLocalizedAdmitCorr ; ./LocalizedAdmitCorr < input.txt
	 @echo
	 @echo SHMAG
	 cd SHMag ; ./SHMag
	 @echo
	 @echo TEST SHROTATE
	 cd SHRotate ; ./TestSHRotate < input.txt
	 @echo
	 @echo "*** The following 8 timing and accuracy tests will take about an hour to complete ***"
	 @echo
	 @echo TIMING ACCURACY GLQ
	 cd TimingAccuracy ; ./TimingAccuracyGLQ < input1 > output1
	 @echo
	 @echo TIMING ACCURACY GLQ
	 cd TimingAccuracy ; ./TimingAccuracyGLQ < input2 > output2
	 @echo
	 @echo TIMING ACCURACY DH
	 cd TimingAccuracy ; ./TimingAccuracyDH < input3 > output3
	 @echo
	 @echo TIMING ACCURACY DH
	 cd TimingAccuracy ; ./TimingAccuracyDH < input4 > output4
	 @echo
	 @echo TIMING ACCURACY DHC
	 cd TimingAccuracy ; ./TimingAccuracyDHC < input5 > output5
	 @echo
	 @echo TIMING ACCURACY DHC
	 cd TimingAccuracy ; ./TimingAccuracyDHC < input6 > output6
	 @echo
	 @echo TIMING ACCURACY GLQC
	 cd TimingAccuracy ; ./TimingAccuracyGLQC < input7 > output7
	 @echo
	 @echo TIMING ACCURACY GLQC
	 cd TimingAccuracy ; ./TimingAccuracyGLQC < input8 > output8

getflags:

ifeq ($(F95),f95)
F95FLAGS ?= -m64 -O3 -YEXT_NAMES=LCS -YEXT_SFX=_ -march=host -speed_math=11
MODFLAG = -p ../modules/
endif

ifeq ($(F95),gfortran)
F95FLAGS ?= -m64 -O3 -march=native -ffast-math
MODFLAG = -I../modules/
endif

ifeq ($(F95),ifort)
F95FLAGS ?= -m64 -free -O3 -Tf
MODFLAG = -I../modules/
endif

ifeq ($(F95),g95)
F95FLAGS ?= -O3 -fno-second-underscore 
MODFLAG = -I../modules/
endif

ifeq ($(F95),pgf90)
F95FLAGS ?= -fast 
MODFLAG = -Imodpath
endif

ifeq ($(origin F95FLAGS), undefined)
F95FLAGS = -m64 -O3
MODFLAG = -I../modules/
endif


clean:
	@rm -f fftw_wisdom/create_wisdom
	@rm -f MarsCrustalThickness/MarsCrustalThickness
	@rm -f MarsCrustalThickness/*.sh
	@rm -f MarsCrustalThickness/*.dat
	@rm -f SHCilmPlus/TestCilmPlus
	@rm -f SHCilmPlus/*.dat
	@rm -f SHExpandDH/TestExpandDH
	@rm -f SHExpandLSQ/TestSHExpandLSQ
	@rm -f SHLocalizedAdmitCorr/LocalizedAdmitCorr
	@rm -f SHLocalizedAdmitCorr/*.dat
	@rm -f SHMag/SHMag
	@rm -f SHMag/*.dat
	@rm -f SHRotate/TestSHRotate
	@rm -f SHRotate/*.sh
	@rm -f SHRotate/*.dat
	@rm -f TimingAccuracy/TimingAccuracyDH
	@rm -f TimingAccuracy/TimingAccuracyDHC
	@rm -f TimingAccuracy/TimingAccuracyGLQ
	@rm -f TimingAccuracy/TimingAccuracyGLQC
	@rm -f TimingAccuracy/output*
	@rm -f TimingAccuracy/*.maxerror
	@rm -f TimingAccuracy/*.rmserror
	@rm -f TimingAccuracy/*.timef
	@rm -f TimingAccuracy/*.timei
	


