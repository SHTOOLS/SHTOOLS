###################################################################################
#
#	make all
#		Compile program in the current directory. Optionally, one
#		can specify the parameters F95="my compiler" and 
#		F95FLAGS="my compiler flags". The default is to use "f95".
#
#	make all2
#		Variant of make all: LAPACK subroutine names have
#		an underscore appended to them in the source files in order to 
#		use FFTW and LAPACK	libraries with conflicting underscore 
#		conventions.
#
#	make all3
#		Variant of make all: FFTW subroutine names have
#		an underscore appended to them in the source files in order to 
#		use FFTW and LAPACK	libraries with conflicting underscore 
#		conventions.
#
#	make python
#		Compile the python wrapper with the f2py compiler. This should be 
#		after the Fortran files are compiled.
#
#	make clean
#		Remove the compiled lib, module, object, and Python files. Also removes 
#		compiled fortran and Python tests.
#
#	make fortran-tests
#		Compile and run example programs and test suite. Optionally, one
#		can specify the parameters F95="my compiler" and 
#		F95FLAGS="my compiler flags", which should be identical to
#		those used to make "all".
#
#	make run-fortran-tests
#		Run all fortran examples and test suite.
#
#	make clean-fortran-tests
#		Delete compiled test suite programs.
#
#   make python-tests
#		Run all python tests
#
#	make clean-python-tests
#		Detelet all compiled python tests
#	make doc
#		Create the man and html-man pages from input POD files.
#		These are PRE-MADE in the distribution, so it shouldn't
#		be necessary to recreate these unless there is some kind
#		of problem.
#
#	make remove-doc
#		Remove the man and html-man pages.
#
#
#	Written by Mark Wieczorek (July 2012).
#
#####################################################################################

VERSION = 3.0

F95 = gfortran
F2PY = f2py
PYTHON = python

FFTW = -lfftw3
LAPACK = -llapack 
BLAS = -lblas

SHELL  = /bin/tcsh
MAKE   = make
FDOCDIR = src/fdoc
PYDOCDIR = src/pydoc
SRCDIR = src
LIBDIR = lib
INCDIR = modules
FEXDIR  = examples/fortran
PEXDIR = examples/python


.PHONY: all all2 all3 python install doc remove-doc clean getflags fortran-tests clean-fortran-tests run-fortran-tests run-python-tests

all: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!
	@echo	
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -Llibpath -lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo
	@echo where modpath and libpath are replaced with their respective paths.
	@echo ---------------------------------------------------------------------------------------------------
	@echo

all2: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all2 F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!
	@echo
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -Llibpath -lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo
	@echo where modpath and libpath are replaced with their respective paths.
	@echo ---------------------------------------------------------------------------------------------------
	@echo

all3: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all3 F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!
	@echo	
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -Llibpath -lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo
	@echo where modpath and libpath are replaced with their respective paths.
	@echo ---------------------------------------------------------------------------------------------------
	@echo 

python: all
	$(F2PY) -I$(INCDIR) -L$(LIBDIR) --f90flags="$(F95FLAGS)" \
	    -c $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95\
	    -lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	$(F2PY) --f90flags="$(F95FLAGS)" -c $(SRCDIR)/PlanetsConstants.f95 -m _constant
	mv _SHTOOLS.so pyshtools/.
	mv _constant.so pyshtools/.
	mkdir -p pyshtools/doc
	./pyshtools/make_docs.py .
	@echo
	@echo MAKE SUCCESSFUL!
	@echo	
	@echo ---------------------------------------------------------------------------------------------------
	@echo import shtools into Python with:
	@echo
	@echo import pyshtools as shtools
	@echo ---------------------------------------------------------------------------------------------------
	@echo 

getflags:
ifeq ($(F95),f95)
# Default Absoft Pro Fortran flags
F95FLAGS ?= -m64 -O3 -YEXT_NAMES=LCS -YEXT_SFX=_ -fpic
#-march=host
MODFLAG = -p modpath
endif

ifeq ($(F95),gfortran)
# Default gfortran flags
#F95FLAGS ?= -m64 -fPIC -O3 # -march=native
F95FLAGS ?= -m64 -fPIC -O3 # -march=native  # -ggdb
MODFLAG = -Imodpath
#LAPACK = #"-framework accelerate" This will compile and run the fortran code, but will not compile the python library.
endif

ifeq ($(F95),ifort)
# Default intel fortran flags
F95FLAGS ?= -m64 -free -O3 -Tf
MODFLAG = -Imodpath
endif

ifeq ($(F95),g95)
# Default g95 flags.
F95FLAGS ?= -O3 -fno-second-underscore 
MODFLAG = -Imodpath
endif

ifeq ($(F95),pgf90)
# Default pgf90 flags
F95FLAGS ?= -fast 
MODFLAG = -Imodpath
endif

ifeq ($(origin F95FLAGS), undefined)
F95FLAGS = -m64 -O3
MODFLAG = -Imodpath
endif


doc: 
	$(MAKE) -C $(FDOCDIR) -f Makefile VERSION=$(VERSION)
	$(MAKE) -C $(PYDOCDIR) -f Makefile VERSION=$(VERSION)
	@echo DOCUMENTATION SUCCESSFULLY CREATED

remove-doc:
	-rm -f man/man1/*.1
	-rm -f www/man/fortran/*.html
	-rm -f www/man/python/*.html
	@echo
	@echo REMOVED MAN AND HTML-MAN FILES

clean:
	-$(MAKE) -C $(SRCDIR) -f Makefile clean
	-rm -f pyshtools/*.so
	-rm -f pyshtools/*.pyc
	-rm -rf pyshtools/doc
	-$(MAKE) -C $(FEXDIR) -f Makefile clean
	-$(MAKE) -C $(PEXDIR) -f Makefile clean
	@echo
	@echo REMOVED LIB, MODULE, OBJECT FILES, FORTRAN TESTS, COMPILED PYTHON FILES AND TESTS

fortran-tests: getflags
	$(MAKE) -C $(FEXDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS)" FFTW=$(FFTW) LAPACK=$(LAPACK) BLAS=$(BLAS)
	@echo
	@echo MAKE OF FORTRAN TEST SUITE SUCCESSFUL
	@echo
	$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests
	@echo
	@echo RAN ALL FORTRAN EXAMPLE AND TESTS

run-fortran-tests:
	$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests
	@echo
	@echo RAN ALL FORTRAN EXAMPLE AND TESTS

clean-fortran-tests:
	$(MAKE) -C $(FEXDIR) -f Makefile clean
	@echo
	@echo REMOVED FORTRAN TEST SUITE EXECUTABLES AND FILES

python-tests:
	$(MAKE) -C $(PEXDIR) -f Makefile all PYTHON=$(PYTHON)
	@echo
	@echo RAN ALL PYTHON TESTS

clean-python-tests:
	$(MAKE) -C $(PEXDIR) -f Makefile clean
	@echo
	@echo REMOVED PYTHON TEST SUITE EXECUTABLES AND FILES
