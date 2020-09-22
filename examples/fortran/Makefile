###############################################################################
#
#   make all
#       Compile example programs and tests in their directories.
#
#   make run-fortran-tests
#       Run the Fortran examples and test suite
#
#   make run-fortran-tests-no-timing
#       Run the Fortran examples and test suite
#
#   make clean
#       Delete the compile executibles.
#
#
#   NOTES
#
#   This Makefile is intended to be invoked by the main Makefile in the SHTOOLS
#   repository. To use this Makefile independently, it is necessary to provide
#   the following arguments:
#
#       F95 - name of the Fortran 95 compiler
#       F95FLAGS - flags to use with the Fortran 95 compiler
#       LIBNAME - name of the SHTOOLS library (without lib and .a)
#       LIBPATH - path to the SHTOOLS library
#       MODFLAG - flag used to identify the path of the compiled modules
#       FFTW - all flags necessary to find and link to the FFTW library
#       LAPACK - all flags necessary to find and link to the LAPACK library
#       BLAS - all flags necessary to find and link to the BLAS library
#
###############################################################################

SHELL = /bin/sh

.PHONY: all clean run-fortran-tests


EXAMPLES = \
	MarsCrustalThickness/MarsCrustalThickness \
	SHCilmPlus/TestCilmPlus \
	SHExpandDH/TestExpandDH \
	SHExpandLSQ/TestSHExpandLSQ \
	SHLocalizedAdmitCorr/LocalizedAdmitCorr \
	SHMag/SHMag \
	SHRotate/TestSHRotate \
	TimingAccuracy/TimingAccuracyDH \
	TimingAccuracy/TimingAccuracyDHC \
	TimingAccuracy/TimingAccuracyGLQ \
	TimingAccuracy/TimingAccuracyGLQC


all: $(EXAMPLES)

%: %.f95 $(LIBPATH)/lib$(LIBNAME).a
	$(F95) $(F95FLAGS) $< $(MODFLAG) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) $(LAPACK) $(BLAS) -o $@

run-fortran-tests: $(EXAMPLES)
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


run-fortran-tests-no-timing: $(EXAMPLES)
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


clean:
	@-rm -f $(EXAMPLES)
	@-rm -f $(addsuffix *.sh, $(dir $(EXAMPLES)))
	@-rm -f $(addsuffix *.dat, $(dir $(EXAMPLES)))
	@-rm -f TimingAccuracy/output*
	@-rm -f TimingAccuracy/*.maxerror
	@-rm -f TimingAccuracy/*.rmserror
	@-rm -f TimingAccuracy/*.timef
	@-rm -f TimingAccuracy/*.timei
