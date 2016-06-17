###################################################################################
#
#	INSTRUCTIONS
#
#	The normal user should only have to use the following commands
#	
#		make, make all		: install both the fortran and python 2 components
#		make fortran		: install the fortan components
#		make fortran-mp		: install the fortan OpenMP components
#		make python			: install the python components
#		make python2		: install the python 2 components
#		make python3		: install the python 3 components
#		make fortran-tests	: compile and run the fortran test/example suite
#		make fortran-tests-mp	: compile and run the fortran test/example suite with OpenMp support
#		make python-tests	: run the python test/example suite
#		make python2-tests	: run the python 2 test/example suite
#		make python3-tests	: run the python 3 test/example suite
#		make install		: place the compiled libraries and docs in /usr/local
#		make uninstall		: remove files copied to /usr/local
#		make clean			: return the folder to its original state
#
#	In some cases, where there are underscore problems when linking to the 
#	LAPACK, BLAS, and FFTW3 libraries, it might be necessary to set
#	LAPACK_UNDERSCORE=1 or FFTW3_UNDERSCORE=1. ALL OF THE OTHER MAKES LISTED
#	BELOW ARE FOR DEVELOPERS ONLY.
#
#	This Makefile accepts the following optional arguments that can be passed
#	using the syntax : make F95="name of f95 compile"
#	
#		F95			: Name and path of the fortran-95 compiler
#		F95FLAGS	: Fortran-95 compiler flags
#		OPENMPFLAGS	: Fortran-95 OpenMP compiler flags
#		F2PY		: Name (including path) of the f2py executable
#		F2PY3		: Name (including path) of the f2py3 executable
#		PYTHON		: Name (including path) of the python executable
#		PYTHON3		: Name (including path) of the python3 executable
#		FFTW		: Name and path of the FFTW3 library of the form "-Lpath -lfftw3"
#		LAPACK		: Name and path of the LAPACK library of the form "-Lpath -llapack"
#		BLAS		: Name and path of the BLAS library of the form "-Lpath -lblas"
#
#
#	LIST OF ALL SUPPORTED MAKE TARGETS
#
#	make, make all
#		Compile fortran and python programs in the current directory.
#
#	make fortran
#		Compile only the fortran 95 component of the archive.
#
#	make fortran-mp
#		Compile only the fortran 95 component of the archive with OpenMP support.
#	
#	make install-fortran
#		Install only fortram 95 component of SHTOOLS.
#
#	make python
#		Compile the python wrapper with the f2py compiler. This should be 
#		after the Fortran files are compiled. The variable PYTHON_VERSION will
#		determine which wrapper to build. If set to "all", then both 2 and 3
#		will be created. If set to either 2 or 3, then only a wrapper for the
#		specified version is built. You can also build for one or the other
#		directly with the python2 or python3 targets below.
#
#	make python2
#		Compile the python 2 wrapper with the f2py compiler. This should be
#		after the Fortran files are compiled.
#
#	make python3
#		Compile the python 3 wrapper with the f2py3 compiler. This should be
#		after the Fortran files are compiled.
#
#	make clean
#		Remove the compiled lib, module, object, and Python files. Also removes 
#		compiled fortran and Python tests.
#
#	make fortran-tests
#		Compile and run example programs and test suite. Optional parameters 
#		should be identical to those used to make "all".
#
#	make run-fortran-tests
#		Run all fortran examples and test suite.
#
#	make clean-fortran-tests
#		Delete compiled test suite programs.
#
#	make python-tests
#		Run all python tests
#
#	make python2-tests
#		Run all python 2 tests
#
#	make python3-tests
#		Run all python 3 tests
#
#	make clean-python-tests
#		Detele all compiled python tests
#
#	make doc
#		Create the man and html-man pages from input Markdown files.
#		These are PRE-MADE in the distribution. To remake these
#		files, it will be necessary to install "pandoc", 
#		"ghc" and "cabal-install" (all using brew on OSX),
#		and then execute "cabal update" and "cabal install pandoc-types".
#
#	make remove-doc
#		Remove the man and html-man pages.
#
#
#	Written by Mark Wieczorek (March 2015).
#
#####################################################################################

VERSION = 3.2
LIBNAME = SHTOOLS
LIBNAMEMP = SHTOOLS-mp

F95 = gfortran
F2PY = f2py
F2PY3 = f2py3
PYTHON = python
PYTHON3 = python3
PYTHON_VERSION = all

SYSLIBPATH = /usr/local/lib

FFTW = -L$(SYSLIBPATH) -lfftw3
FFTW_UNDERSCORE = 0
LAPACK = -llapack 
LAPACK_UNDERSCORE = 0
BLAS = -lblas

SHELL = /bin/sh
FDOCDIR = src/fdoc
PYDOCDIR = src/pydoc
SRCDIR = src
LIBDIR = lib
INCDIR = modules
FEXDIR = examples/fortran
PEXDIR = examples/python

LIBPATH = $(PWD)/$(LIBDIR)
MODPATH = $(PWD)/$(INCDIR)
PYPATH = $(PWD)
SYSMODPATH = /usr/local/include
SYSPYPATH = $(shell $(PYTHON) -c 'import sysconfig; print(sysconfig.get_path("platlib"))')
SYSPY3PATH := $(shell $(PYTHON3) -c 'import sysconfig; print(sysconfig.get_path("platlib"))')
PY3EXT := $(shell $(PYTHON3) -c 'import sysconfig; print(sysconfig.get_config_var("SO"))')
SYSSHAREPATH =/usr/local/share
SYSDOCPATH = /usr/local/share/doc

ifeq ($(F95),f95)
# Default Absoft f95 flags
F95FLAGS ?= -m64 -O3 -YEXT_NAMES=LCS -YEXT_SFX=_ -fpic -speed_math=10
#-march=host
MODFLAG = -p $(MODPATH)
SYSMODFLAG = -p $(SYSMODPATH)
OPENMPFLAGS ?= -openmp
else ifeq ($(F95),gfortran)
# Default gfortran flags
F95FLAGS ?= -m64 -fPIC -O3 -ffast-math
# -march=native
MODFLAG = -I$(MODPATH)
SYSMODFLAG = -I$(SYSMODPATH)
OPENMPFLAGS ?= -fopenmp
else ifeq ($(F95),ifort)
# Default intel fortran flags
F95FLAGS ?= -m64 -free -O3 -Tf
MODFLAG = -I$(MODPATH)
SYSMODFLAG = -I$(SYSMODPATH)
OPENMPFLAGS ?=
else ifeq ($(F95),g95)
# Default g95 flags.
F95FLAGS ?= -O3 -fno-second-underscore
MODFLAG = -I$(MODPATH)
SYSMODFLAG = -I$(SYSMODPATH)
OPENMPFLAGS ?=
else ifeq ($(F95),pgf90)
# Default pgf90 flags
F95FLAGS ?= -fast 
MODFLAG = -Imodpath
SYSMODFLAG = -Imodpath
OPENMPFLAGS ?=
else ifeq ($(origin F95FLAGS), undefined)
F95FLAGS = -m64 -O3
MODFLAG = -I$(MODPATH)
SYSMODFLAG = -I$(SYSMODPATH)
OPENMPFLAGS ?=
endif

ifeq ($(LAPACK_UNDERSCORE),1)
LAPACK_FLAGS = -DLAPACK_UNDERSCORE
endif

ifeq ($(FFTW3_UNDERSCORE),1)
FFTW3_FLAGS = -DFFTW3_UNDERSCORE
endif


.PHONY: all fortran python python2 python3 install doc remove-doc clean\
	fortran-tests clean-fortran-tests run-fortran-tests run-python-tests run-python2-tests run-python3-tests\
	install-fortran install-python install-python2 install-python3 uninstall fortran-mp


all: fortran python

fortran:
	$(MAKE) -C $(SRCDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS)" LIBNAME="$(LIBNAME)" LAPACK_FLAGS="$(LAPACK_FLAGS)" FFTW3_FLAGS="$(FFTW3_FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!
	@echo
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your Fortran code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo ---------------------------------------------------------------------------------------------------
	@echo

fortran-mp:
	$(MAKE) -C $(SRCDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS) $(OPENMPFLAGS)" LIBNAME="$(LIBNAMEMP)" LAPACK_FLAGS="$(LAPACK_FLAGS)" FFTW3_FLAGS="$(FFTW3_FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!
	@echo
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your Fortran code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) $(OPENMPFLAGS) -L$(LIBPATH) -l$(LIBNAMEMP) $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo ---------------------------------------------------------------------------------------------------
	@echo

ifeq ($(PYTHON_VERSION),all)
python: python2 python3
install-python: install-python2 install-python3
python-tests: python2-tests python3-tests
clean-python-tests: clean-python2-tests clean-python3-tests
else ifeq ($(PYTHON_VERSION),2)
python: python2
install-python: install-python2
python-tests: python2-tests
clean-python-tests: clean-python2-tests
else ifeq ($(PYTHON_VERSION),3)
python: python3
install-python: install-python3
python-tests: python3-tests
clean-python-tests: clean-python3-tests
else
$(error $(PYTHON_VERSION) is unsupported.)
endif

python2: pyshtools/_SHTOOLS.so pyshtools/_constant.so
	mkdir -p pyshtools/doc
	./pyshtools/make_docs.py .
	@echo
	@echo MAKE SUCCESSFUL!
	@echo
	@echo ---------------------------------------------------------------------------------------------------
	@echo import shtools into Python 2 with:
	@echo
	@echo import sys
	@echo sys.path.append\(\'$(PYPATH)\'\)
	@echo import pyshtools as shtools
	@echo ---------------------------------------------------------------------------------------------------
	@echo

python3: pyshtools/_SHTOOLS$(PY3EXT) pyshtools/_constant$(PY3EXT)
	mkdir -p pyshtools/doc
	$(PYTHON3) ./pyshtools/make_docs.py .
	@echo
	@echo MAKE SUCCESSFUL!
	@echo
	@echo ---------------------------------------------------------------------------------------------------
	@echo import shtools into Python 3 with:
	@echo
	@echo import sys
	@echo sys.path.append\(\'$(PYPATH)\'\)
	@echo import pyshtools as shtools
	@echo ---------------------------------------------------------------------------------------------------
	@echo

pyshtools/_SHTOOLS.so: $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95 $(LIBDIR)/lib$(LIBNAME).a
	$(F2PY) --quiet -I$(INCDIR) -L$(LIBDIR) --f90flags="$(F95FLAGS)" \
		-c $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95\
		-l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS)
	mv _SHTOOLS.so pyshtools/.

pyshtools/_constant.so: $(SRCDIR)/PlanetsConstants.f95
	$(F2PY) --quiet --f90flags="$(F95FLAGS)" -c $(SRCDIR)/PlanetsConstants.f95 -m _constant
	mv _constant.so pyshtools/.

pyshtools/_SHTOOLS$(PY3EXT): $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95 $(LIBDIR)/lib$(LIBNAME).a
	$(F2PY3) --quiet -I$(INCDIR) -L$(LIBDIR) --f90flags="$(F95FLAGS)" \
		-c $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95\
		-l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS)
	mv _SHTOOLS$(PY3EXT) pyshtools/.

pyshtools/_constant$(PY3EXT): $(SRCDIR)/PlanetsConstants.f95
	$(F2PY3) --quiet --f90flags="$(F95FLAGS)" -c $(SRCDIR)/PlanetsConstants.f95 -m _constant
	mv _constant$(PY3EXT) pyshtools/.

install: install-fortran install-python

install-python2: python2
	mkdir -pv $(SYSPYPATH)/pyshtools
	cp -R $(filter-out %$(PY3EXT), $(wildcard pyshtools/*)) $(SYSPYPATH)/pyshtools/
	@echo ---------------------------------------------------------------------------------------------------
	@echo import shtools into Python 2 with:
	@echo
	@echo import sys
	@echo sys.path.append\(\'$(SYSPYPATH)\'\)
	@echo import pyshtools as shtools
	@echo ---------------------------------------------------------------------------------------------------

install-python3: python3
	mkdir -pv $(SYSPY3PATH)/pyshtools
	cp -R $(filter-out %.so, $(wildcard pyshtools/*)) $(SYSPY3PATH)/pyshtools/
	cp -R $(filter %$(PY3EXT), $(wildcard pyshtools/*)) $(SYSPY3PATH)/pyshtools/
	@echo ---------------------------------------------------------------------------------------------------
	@echo import shtools into Python 3 with:
	@echo
	@echo import sys
	@echo sys.path.append\(\'$(SYSPY3PATH)\'\)
	@echo import pyshtools as shtools
	@echo ---------------------------------------------------------------------------------------------------

uninstall:
	-rm -r $(SYSPYPATH)/pyshtools/
	-rm -r $(SYSPY3PATH)/pyshtools/
	-rm -r $(SYSLIBPATH)/lib$(LIBNAME).a
	-rm -r $(SYSMODPATH)/fftw3.mod
	-rm -r $(SYSMODPATH)/planetsconstants.mod
	-rm -r $(SYSMODPATH)/shtools.mod
	-rm -r $(SYSSHAREPATH)/shtools/examples/
	-rm -r $(SYSDOCPATH)/shtools/index.html
	-rm -r $(SYSDOCPATH)/shtools/www/
	$(MAKE) -C $(FDOCDIR) -f Makefile VERSION=$(VERSION) uninstall
	$(MAKE) -C $(PYDOCDIR) -f Makefile VERSION=$(VERSION) uninstall

install-fortran: fortran
	-rm -r $(SYSMODPATH)/fftw3.mod
	-rm -r $(SYSMODPATH)/planetsconstants.mod
	-rm -r $(SYSMODPATH)/shtools.mod
	mkdir -pv $(SYSLIBPATH)
	cp $(LIBDIR)/lib$(LIBNAME).a $(SYSLIBPATH)/lib$(LIBNAME).a
	mkdir -pv $(SYSMODPATH)
	cp $(INCDIR)/*.mod $(SYSMODPATH)/
	mkdir -pv $(SYSSHAREPATH)/shtools
	cp -R examples $(SYSSHAREPATH)/shtools/
	mkdir -pv $(SYSSHAREPATH)/man/man1
	cp -R man/man1/ $(SYSSHAREPATH)/man/man1/
	mkdir -pv $(SYSDOCPATH)/shtools
	cp index.html $(SYSDOCPATH)/shtools/index.html
	cp -R www $(SYSDOCPATH)/shtools/
	awk '{gsub("../../lib","/usr/local/lib");print}' "examples/Makefile" > "temp.txt"
	awk '{gsub("../../modules","/usr/local/include");print}' "temp.txt" > "temp2.txt"
	cp temp2.txt "/usr/local/share/shtools/examples/Makefile"
	awk '{gsub("../../lib","/usr/local/lib");print}' "examples/fortran/Makefile" > "temp.txt"
	awk '{gsub("../../modules","/usr/local/include");print}' "temp.txt" > "temp2.txt"
	cp temp2.txt "$(SYSSHAREPATH)/shtools/examples/fortran/Makefile"
	rm temp.txt
	rm temp2.txt
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your Fortran code with the following flags:
	@echo
	@echo $(F95) $(SYSMODFLAG) $(F95FLAGS) -L$(SYSLIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo ---------------------------------------------------------------------------------------------------
	@echo 

doc: 
	@$(MAKE) -C $(FDOCDIR) -f Makefile VERSION=$(VERSION)
	@$(MAKE) -C $(PYDOCDIR) -f Makefile VERSION=$(VERSION)
	@echo DOCUMENTATION SUCCESSFULLY CREATED

remove-doc:
	-rm -f man/man1/*.1
	-rm -f www/man/fortran/*.html
	-rm -f www/man/python/*.html
	@echo
	@echo REMOVED MAN AND HTML-MAN FILES

clean: clean-libs clean-fortran-tests clean-python-tests

clean-libs:
	-$(MAKE) -C $(SRCDIR) -f Makefile clean
	-rm -f lib/lib$(LIBNAME).a
	-rm -f lib/lib$(LIBNAMEMP).a
	-rm -f pyshtools/*.so
	-rm -f pyshtools/*.pyc
	-rm -rf pyshtools/__pycache__/
	-rm -rf pyshtools/doc
	@echo
	@echo REMOVED LIB, MODULE, OBJECT FILES, COMPILED PYTHON FILES AND TESTS

fortran-tests: fortran
	$(MAKE) -C $(FEXDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS)" LIBNAME="$(LIBNAME)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo MAKE OF FORTRAN TEST SUITE SUCCESSFUL
	@echo
	$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests
	@echo
	@echo RAN ALL FORTRAN EXAMPLE AND TESTS

fortran-tests-mp: fortran-mp
	$(MAKE) -C $(FEXDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS) $(OPENMPFLAGS)" LIBNAME="$(LIBNAMEMP)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo MAKE OF FORTRAN TEST SUITE SUCCESSFUL
	@echo
	$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests
	@echo
	@echo RAN ALL FORTRAN EXAMPLE AND TESTS

run-fortran-tests: fortran
	$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests
	@echo
	@echo RAN ALL FORTRAN EXAMPLE AND TESTS

clean-fortran-tests:
	$(MAKE) -C $(FEXDIR) -f Makefile clean
	@echo
	@echo REMOVED FORTRAN TEST SUITE EXECUTABLES AND FILES

python2-tests: python2
	$(MAKE) -C $(PEXDIR) -f Makefile all PYTHON=$(PYTHON)
	@echo
	@echo RAN ALL PYTHON TESTS

python3-tests: python3
	$(MAKE) -C $(PEXDIR) -f Makefile all PYTHON=$(PYTHON3)
	@echo
	@echo RAN ALL PYTHON 3 TESTS

clean-python2-tests:
	$(MAKE) -C $(PEXDIR) -f Makefile clean
	@echo
	@echo REMOVED PYTHON TEST SUITE EXECUTABLES AND FILES

clean-python3-tests:
	$(MAKE) -C $(PEXDIR) -f Makefile clean
	@echo
	@echo REMOVED PYTHON TEST SUITE EXECUTABLES AND FILES
