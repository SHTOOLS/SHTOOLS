###################################################################################
#
#	INSTRUCTIONS
#
#	The normal user should only have to use the following commands
#	
#		make, make all		: install both the fortran and python components
#		make fortran		: install the fortan components
#		make python			: install the python components	   
#		make fortran-tests	: compile and run the fortran test/example suite
#		make python-tests	: run the python test/example suite
#		make install		: place the compiled libraries and docs in /usr/local
#		make uninstall		: remove files copied to /usr/local
#		make clean			: return the folder to its original state
#
#	In some cases, where there are underscore problems when linking to the 
#	LAPACK, BLAS, and FFTW3 libraries, it might be necessary to use
#	make all2 or make all3 instead of make all. ALL OF THE OTHER
#	MAKES LISTED BELOW ARE FOR DEVELOPERS ONLY.
#
#	This Makefile accepts the following optional arguments that can be passed
#	using the syntax : make F95="name of f95 compile"
#	
#		F95			: Name and path of the fortran-95 compiler
#		F95FLAGS	: Fortran-95 compiler flags
#		F2PY		: Name (including path) of the f2py executable
#		PYTHON		: Name (including path) of the python executable
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
#	make fortran2
#		Variant of make all: LAPACK subroutine names have
#		an underscore appended to them in the source files in order to 
#		use FFTW and LAPACK libraries with conflicting underscore 
#		conventions.
#
#	make fortran3
#		Variant of make all: FFTW subroutine names have
#		an underscore appended to them in the source files in order to 
#		use FFTW and LAPACK libraries with conflicting underscore 
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
#	make clean-python-tests
#		Detele all compiled python tests
#
#	make doc
#		Create the man and html-man pages from input POD files.
#		These are PRE-MADE in the distribution. To remake these
#		files, it will be necessary to install "pandoc", 
#       "ghc" and "cabal-install" (all using brew on OSX),
#		and then execute "cabal update| and "cabal install pandoc-types".
#
#	make remove-doc
#		Remove the man and html-man pages.
#
#
#	Written by Mark Wieczorek (March 2015).
#
#####################################################################################

VERSION = 3.0

F95 = gfortran
F2PY = f2py
PYTHON = python

FFTW = -lfftw3
LAPACK = -llapack 
BLAS = -lblas

SHELL = /bin/tcsh
MAKE = make
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
SYSLIBPATH = /usr/local/lib
SYSPYPATH = /usr/local/lib/python2.7/site-packages
SYSSHAREPATH =/usr/local/share
SYSDOCPATH = /usr/local/share/doc

.PHONY: all all2 all3 fortran fortran2 fortran3 python install doc remove-doc clean\
	getflags fortran-tests clean-fortran-tests run-fortran-tests run-python-tests\
	install uninstall

all: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE OF FORTRAN COMPONENT SUCCESSFUL!
	@echo
	$(F2PY) --quiet -I$(INCDIR) -L$(LIBDIR) --f90flags="$(F95FLAGS)" \
		-c $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95\
		-lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	$(F2PY) --quiet --f90flags="$(F95FLAGS)" -c $(SRCDIR)/PlanetsConstants.f95 -m _constant
	mv _SHTOOLS.so pyshtools/.
	mv _constant.so pyshtools/.
	mkdir -p pyshtools/doc
	./pyshtools/make_docs.py .
	@echo
	@echo MAKE OF PYTHON COMPONENT SUCCESSFUL!
	@echo 
	@echo -----------------------------------------------------------------------------------------
	@echo Compile your Fortran code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo
	@echo import shtools into Python with:
	@echo
	@echo import sys
	@echo sys.path.append\(\'$(PYPATH)\'\)
	@echo import pyshtools as shtools
	@echo --------------------------------------------------------------------------------------------
	@echo

all2: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all2 F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE OF FORTRAN COMPONENT SUCCESSFUL!
	@echo
	$(F2PY) --quiet -I$(INCDIR) -L$(LIBDIR) --f90flags="$(F95FLAGS)" \
		-c $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95\
		-lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	$(F2PY) --quiet --f90flags="$(F95FLAGS)" -c $(SRCDIR)/PlanetsConstants.f95 -m _constant
	mv _SHTOOLS.so pyshtools/.
	mv _constant.so pyshtools/.
	mkdir -p pyshtools/doc
	./pyshtools/make_docs.py .
	@echo
	@echo MAKE OF PYTHON COMPONENT SUCCESSFUL!
	@echo
	@echo ------------------------------------------------------------------------------------------
	@echo Compile your Fortran code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo
	@echo import shtools into Python with:
	@echo
	@echo import sys
	@echo sys.path.append\(\'$(PYPATH)\'\)
	@echo import pyshtools as shtools
	@echo -------------------------------------------------------------------------------------------
	@echo	
	
all3: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all3 F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE OF FORTRAN COMPONENT SUCCESSFUL!
	@echo 
	$(F2PY) --quiet -I$(INCDIR) -L$(LIBDIR) --f90flags="$(F95FLAGS)" \
		-c $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95\
		-lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	$(F2PY) --quiet --f90flags="$(F95FLAGS)" -c $(SRCDIR)/PlanetsConstants.f95 -m _constant
	mv _SHTOOLS.so pyshtools/.
	mv _constant.so pyshtools/.
	mkdir -p pyshtools/doc
	./pyshtools/make_docs.py .
	@echo
	@echo MAKE OF PYTHON COMPONENT SUCCESSFUL!
	@echo 
	@echo ---------------------------------------------------------------------------------------------
	@echo Compile your Fortran code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo
	@echo import shtools into Python with:
	@echo
	@echo import sys
	@echo sys.path.append\(\'$(PYPATH)\'\)
	@echo import pyshtools as shtools
	@echo ----------------------------------------------------------------------------------------------
	@echo

fortran: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!
	@echo
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your Fortran code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo
	@echo ---------------------------------------------------------------------------------------------------
	@echo
	
fortran2: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all2 F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!
	@echo
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your Fortran code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo
	@echo ---------------------------------------------------------------------------------------------------
	@echo

fortran3: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all3 F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!
	@echo
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your Fortran code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo
	@echo ---------------------------------------------------------------------------------------------------
	@echo 

python: getflags
	$(F2PY) --quiet -I$(INCDIR) -L$(LIBDIR) --f90flags="$(F95FLAGS)" \
		-c $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95\
		-lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	$(F2PY) --quiet --f90flags="$(F95FLAGS)" -c $(SRCDIR)/PlanetsConstants.f95 -m _constant
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
	@echo import sys
	@echo sys.path.append\(\'$(PYPATH)\'\)
	@echo import pyshtools as shtools
	@echo ---------------------------------------------------------------------------------------------------
	@echo 

install: getflags
	mkdir -p $(SYSPYPATH)
	cp -r pyshtools/ $(SYSPYPATH)/pyshtools/
	mkdir -p $(SYSLIBPATH)
	cp $(LIBDIR)/libSHTOOLS.a $(SYSLIBPATH)/libSHTOOLS.a
	mkdir -p $(SYSMODPATH)
	cp $(INCDIR)/fftw3.mod $(SYSMODPATH)/fftw3.mod
	cp $(INCDIR)/planetsconstants.mod $(SYSMODPATH)/planetsconstants.mod
	cp $(INCDIR)/shtools.mod $(SYSMODPATH)/shtools.mod
	mkdir -p $(SYSSHAREPATH)/shtools
	cp -r examples/ $(SYSSHAREPATH)/shtools/examples/
	mkdir -p $(SYSSHAREPATH)/man/man1
	cp -r man/man1/ $(SYSSHAREPATH)/man/man1/
	mkdir -p $(SYSDOCPATH)/shtools
	cp index.html $(SYSDOCPATH)/shtools/index.html
	cp -r www/ $(SYSDOCPATH)/shtools/www/
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
	@echo $(F95) $(SYSMODFLAG) $(F95FLAGS) -L$(SYSLIBPATH) -lSHTOOLS $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo
	@echo import shtools into Python with:
	@echo
	@echo import sys
	@echo sys.path.append\(\'$(SYSPYPATH)\'\)
	@echo import pyshtools as shtools
	@echo ---------------------------------------------------------------------------------------------------

uninstall:
	-rm -r $(SYSPYPATH)/pyshtools/
	-rm -r $(SYSLIBPATH)/libSHTOOLS.a
	-rm -r $(SYSMODPATH)/fftw3.mod
	-rm -r $(SYSMODPATH)/planetsconstants.mod
	-rm -r $(SYSMODPATH)/shtools.mod
	-rm -r $(SYSSHAREPATH)/shtools/examples/
	-rm -r $(SYSDOCPATH)/shtools/index.html
	-rm -r $(SYSDOCPATH)/shtools/www/
	$(MAKE) -C $(FDOCDIR) -f Makefile VERSION=$(VERSION) uninstall
	$(MAKE) -C $(PYDOCDIR) -f Makefile VERSION=$(VERSION) uninstall


getflags:
ifeq ($(F95),f95)
# Default Absoft f95 flags
F95FLAGS ?= -m64 -O3 -YEXT_NAMES=LCS -YEXT_SFX=_ -fpic -speed_math=11
#-march=host
MODFLAG = -p $(MODPATH)
SYSMODFLAG = -p $(SYSMODPATH)
endif

ifeq ($(F95),gfortran)
# Default gfortran flags
F95FLAGS ?= -m64 -fPIC -O3 -ffast-math
# -march=native
MODFLAG = -I$(MODPATH)
SYSMODFLAG = -I$(SYSMODPATH)
endif

ifeq ($(F95),ifort)
# Default intel fortran flags
F95FLAGS ?= -m64 -free -O3 -Tf
MODFLAG = -I$(MODPATH)
SYSMODFLAG = -I$(SYSMODPATH)
endif

ifeq ($(F95),g95)
# Default g95 flags.
F95FLAGS ?= -O3 -fno-second-underscore 
MODFLAG = -I$(MODPATH)
SYSMODFLAG = -I$(SYSMODPATH)
endif

ifeq ($(F95),pgf90)
# Default pgf90 flags
F95FLAGS ?= -fast 
MODFLAG = -Imodpath
SYSMODFLAG = -Imodpath
endif

ifeq ($(origin F95FLAGS), undefined)
F95FLAGS = -m64 -O3
MODFLAG = -I$(MODPATH)
SYSMODFLAG = -I$(SYSMODPATH)
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

