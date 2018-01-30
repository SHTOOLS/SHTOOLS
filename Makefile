###############################################################################
#
#   INSTRUCTIONS
#
#   The normal user should only have to use the following commands
#
#       make (all, fortran) : Install the Fortran 95 components.
#       make fortran-mp     : Install the Fortan OpenMP components.
#       make python         : Install the Python components (versions
#                             determined by PYTHON_VERSION).
#       make python2        : Install the Python 2 components.
#       make python3        : Install the Python 3 components.
#       make fortran-tests  : Compile and run the Fortran test/example suite.
#       make fortran-tests-mp : Compile and run the fortran test/example suite
#                               with OpenMp support.
#       make fortran-tests-no-timing    : Compile and run the Fortran
#                                         test/example suite (excluding timing
#                                         tests).
#       make fortran-tests-no-timing-mp : Compile and run the Fortran OpenMP
#                                         test/example suite (excluding timing
#                                         tests).
#       make python-tests   : Run the Python test/example suite (versions
#                             determined by PYTHON_VERSION).
#       make python2-tests  : Run the Python 2 test/example suite.
#       make python3-tests  : Run the Python 3 test/example suite.
#       make install        : Place the compiled libraries and docs in
#                             $(DESTDIR)$(PREFIX) [default is /usr/local].
#       make uninstall      : Remove files copied to $(DESTDIR)$(PREFIX).
#       make clean          : Return the folder to its original state.
#
#
#   In some cases, where there are underscore problems when linking to the
#   LAPACK, BLAS, and FFTW3 libraries, it might be necessary to set one of
#   the following:
#
#       LAPACK_UNDERSCORE = 1
#       FFTW3_UNDERSCORE = 1
#
#
#   This Makefile accepts the following optional arguments that can be passed
#   using the syntax : make NAME="options"
#
#       F95         : Name and path of the Fortran-95 compiler.
#       F95FLAGS    : Fortran-95 compiler flags.
#       OPENMPFLAGS : Fortran-95 OpenMP compiler flags.
#       F2PY        : Name (including path) of the Python 2 f2py executable.
#       F2PY3       : Name (including path) of the Python 3 f2py executable.
#       PYTHON      : Name (including path) of the Python 2 executable.
#       PYTHON3     : Name (including path) of the Python 3 executable.
#       PYTHON_VERSION : Which Python versions of shtools to install: either
#                        2, 3, or all.
#       FFTW        : Name and path of the FFTW3 library of the form
#                     "-Lpath -lfftw3".
#       LAPACK      : Name and path of the LAPACK library of the form
#                     "-Lpath -llapack"
#       BLAS        : Name and path of the BLAS library of the form
#                     "-Lpath -lblas"
#       PREFIX      : Prefix appended to the system lib, share, doc, and
#                     include directories [default /usr/local].
#       DESTDIR     : Path appended to all system install directories for use
#                     with staged installs [default is none].
#
#
#   LIST OF ALL SUPPORTED MAKE TARGETS
#
#   make, make all, make fortran
#       Compile the Fortran 95 components the current directory.
#
#   make fortran-mp
#       Compile only the fortran 95 component of the archive with OpenMP
#       support.
#
#   make install-fortran
#       Install only the fortran 95 component of SHTOOLS.
#
#   make python
#       Compile the Python components of SHTOOLS with the f2py compiler. This
#       should be done after the Fortran files are compiled. The variable
#       PYTHON_VERSION will determine which versions to build. If set to "all",
#       then both Python versions 2 and 3 will be created. If set to either 2
#       or 3, then only the specified version is built. You can also build for
#       one or the other directly with the python2 or python3 targets below.
#
#   make python2
#       Compile the Python 2 SHTOOLS components with the f2py compiler.
#
#   make python3
#       Compile the python 3 SHTOOLS components with the f2py3 compiler.
#
#   make clean
#       Remove the compiled lib, module, object, and Python files. Also removes
#      compiled fortran and Python tests.
#
#   make fortran-tests
#       Compile and run example Fortran programs and test suite. Optional
#       parameters should be identical to those used to make "all".
#
#   make run-fortran-tests
#       Run all Fortran examples and test suite.
#
#   make run-fortran-tests-no-timing
#       Run all Fortran examples and test suite (excluding the timing tests).
#
#   make clean-fortran-tests
#       Delete compiled Fortran test suite programs.
#
#   make python-tests
#       Run all Python tests, where the variable PYTHON_VERSION determines
#       which versions to run.
#
#   make python2-tests
#       Run all Python 2 tests.
#
#   make python3-tests
#       Run all Python 3 tests.
#
#   make python-tests-no-timing
#       Run all Python tests (with exception of timing/accuracy tests), where
#       the variable PYTHON_VERSION determines which versions to run.
#
#   make python2-tests-no-timing
#       Run all Python tests (with exception of timing/accuracy tests).
#
#   make python3-tests-no-timing
#       Run all Python tests (with exception of timing/accuracy tests).
#
#   make clean-python-tests
#       Detele all compiled Python test files.
#
#   make notebooks
#       Run notebooks and convert to html for web documentation.
#
#   make remove-notebooks
#       Remove html notebooks.
#
#   make doc
#       Create the man pages from input markdown files and create the static
#       web site. Both of these are PRE-MADE in the distribution. To remake
#       these files, it will be necessary to install "pandoc", "ghc" and
#       "cabal-install" (all using brew on OSX), and then execute
#       "cabal update" and "cabal install pandoc-types".
#
#   make remove-doc
#       Remove the man and html-man pages.
#
#   make www
#       Make the static html web documention in the directory www using Jekyll.
#
#   make remove-www
#       Remove the directory containing the static html web site.
#
###############################################################################

VERSION = 4.2
LIBNAME = SHTOOLS
LIBNAMEMP = SHTOOLS-mp

F95 = gfortran
F2PY = f2py
F2PY3 = f2py3
PYTHON = python
PYTHON3 = python3
PYTHON_VERSION = all
JUPYTER = "jupyter nbconvert --ExecutePreprocessor.kernel_name=python2"
JUPYTER3 = "jupyter nbconvert --ExecutePreprocessor.kernel_name=python3"
JEKYLL = bundle exec jekyll

PREFIX = /usr/local
SYSLIBPATH = $(PREFIX)/lib

FFTW = -L$(SYSLIBPATH) -lfftw3
FFTW_UNDERSCORE = 0
LAPACK_UNDERSCORE = 0

SHELL = /bin/sh
FDOCDIR = src/fdoc
PYDOCDIR = src/pydoc
SRCDIR = src
LIBDIR = lib
INCDIR = modules
FEXDIR = examples/fortran
PEXDIR = examples/python
NBDIR = examples/notebooks
WWWSRC = docs
WWWDEST = www

LIBPATH = $(PWD)/$(LIBDIR)
MODPATH = $(PWD)/$(INCDIR)
PYPATH = $(PWD)
SYSMODPATH = $(PREFIX)/include
SYSPYPATH = $(shell $(PYTHON) -c 'import sysconfig; print(sysconfig.get_path("platlib"))')
SYSPY3PATH = $(shell $(PYTHON3) -c 'import sysconfig; print(sysconfig.get_path("platlib"))')
PY3EXT = $(shell $(PYTHON3) -c 'import sysconfig; print(sysconfig.get_config_var("SO"))' || echo nopy3)
SYSSHAREPATH = $(PREFIX)/share
SYSDOCPATH = $(PREFIX)/share/doc

ifeq ($(F95), f95)
# Default Absoft f95 flags
F95FLAGS ?= -m64 -O3 -YEXT_NAMES=LCS -YEXT_SFX=_ -fpic -speed_math=10
#-march=host
MODFLAG = -p $(MODPATH)
SYSMODFLAG = -p $(SYSMODPATH)
OPENMPFLAGS ?= -openmp
BLAS ?= -lblas
LAPACK ?= -llapack
else ifeq ($(F95), gfortran)
# Default gfortran flags
F95FLAGS ?= -m64 -fPIC -O3 -ffast-math
# -march=native
MODFLAG = -I$(MODPATH)
SYSMODFLAG = -I$(SYSMODPATH)
OPENMPFLAGS ?= -fopenmp
BLAS ?= -lblas
LAPACK ?= -llapack
else ifeq ($(F95), ifort)
# Default intel fortran flags
F95FLAGS ?= -m64 -fpp -free -O3 -Tf
MODFLAG = -I$(MODPATH)
SYSMODFLAG = -I$(SYSMODPATH)
OPENMPFLAGS ?= -qopenmp
LAPACK ?= -mkl
BLAS ?= 
else ifeq ($(F95), g95)
# Default g95 flags.
F95FLAGS ?= -O3 -fno-second-underscore
MODFLAG = -I$(MODPATH)
SYSMODFLAG = -I$(SYSMODPATH)
OPENMPFLAGS ?=
BLAS ?= -lblas
LAPACK ?= -llapack
else ifeq ($(F95), pgf90)
# Default pgf90 flags
F95FLAGS ?= -fast 
MODFLAG = -Imodpath
SYSMODFLAG = -Imodpath
OPENMPFLAGS ?=
BLAS ?= -lblas
LAPACK ?= -llapack
else ifeq ($(origin F95FLAGS), undefined)
F95FLAGS = -m64 -O3
MODFLAG = -I$(MODPATH)
SYSMODFLAG = -I$(SYSMODPATH)
OPENMPFLAGS ?=
BLAS ?= -lblas
LAPACK ?= -llapack
endif

ifeq ($(LAPACK_UNDERSCORE),1)
LAPACK_FLAGS = -DLAPACK_UNDERSCORE
endif

ifeq ($(FFTW3_UNDERSCORE),1)
FFTW3_FLAGS = -DFFTW3_UNDERSCORE
endif


.PHONY: all fortran python python2 python3 install doc remove-doc\
	fortran-tests run-fortran-tests python-tests\
	python2-tests python3-tests python-tests-no-timing python2-tests-no-timing\
	python3-tests-no-timing install-fortran install-python\
	install-python2 install-python3 uninstall fortran-mp clean\
	clean-fortran-tests clean-python-tests clean-python2 clean-python3\
	clean-libs remove-notebooks notebooks notebooks2 notebooks3 www remove-www


all: fortran

fortran:
	mkdir -pv lib
	mkdir -pv modules
	$(MAKE) -C $(SRCDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS)" LIBNAME="$(LIBNAME)" LAPACK_FLAGS="$(LAPACK_FLAGS)" FFTW3_FLAGS="$(FFTW3_FLAGS)"
	@echo "--> make fortran successful!"
	@echo
	@echo "Compile your Fortran code with the following flags:"
	@echo "---------------------------------------------------"
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo

fortran-mp:
# Delete .o files before and after compiling with OpenMP to avoid issues with "fortran" build.
	mkdir -pv lib
	mkdir -pv modules
	-$(MAKE) -C $(SRCDIR) -f Makefile clean-obs-mod
	$(MAKE) -C $(SRCDIR) -f Makefile all F95=$(F95) F95FLAGS="$(OPENMPFLAGS) $(F95FLAGS)" LIBNAME="$(LIBNAMEMP)" LAPACK_FLAGS="$(LAPACK_FLAGS)" FFTW3_FLAGS="$(FFTW3_FLAGS)"
	-$(MAKE) -C $(SRCDIR) -f Makefile clean-obs-mod
	@echo "--> make fortran-mp successful!"
	@echo
	@echo "Compile your Fortran code with the following flags:"
	@echo "---------------------------------------------------"
	@echo $(F95) $(MODFLAG) $(OPENMPFLAGS) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAMEMP) $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo

ifeq ($(PYTHON_VERSION),all)
python: python2 python3
install-python: install-python2 install-python3
python-tests: python2-tests python3-tests
python-tests-no-timing: python2-tests-no-timing python3-tests-no-timing
notebooks: notebooks3 notebooks2
else ifeq ($(PYTHON_VERSION),2)
python: python2
install-python: install-python2
python-tests: python2-tests
python-tests-no-timing: python2-tests-no-timing
notebooks: notebooks2
else ifeq ($(PYTHON_VERSION),3)
python: python3
install-python: install-python3
python-tests: python3-tests
python-tests-no-timing: python3-tests-no-timing
notebooks: notebooks3
else
$(error $(PYTHON_VERSION) is unsupported.)
endif

python2: fortran pyshtools/_SHTOOLS.so pyshtools/_constant.so
	mkdir -p pyshtools/doc
	$(PYTHON) ./pyshtools/make_docs.py . .
	@echo
	@echo "--> make python2 successful!"
	@echo
	@echo "Import shtools into Python 2 with:"
	@echo "----------------------------------"
	@echo "import sys"
	@echo "sys.path.append('$(PYPATH)')"
	@echo "import pyshtools"
	@echo

python3: fortran pyshtools/_SHTOOLS$(PY3EXT) pyshtools/_constant$(PY3EXT)
	mkdir -p pyshtools/doc
	$(PYTHON3) ./pyshtools/make_docs.py . .
	@echo
	@echo "--> make python3 successful!"
	@echo
	@echo "Import shtools into Python 3 with:"
	@echo "----------------------------------"
	@echo "import sys"
	@echo "sys.path.append('$(PYPATH)')"
	@echo "import pyshtools"
	@echo

pyshtools/_SHTOOLS.so: $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95 fortran
	$(F2PY) --quiet -I$(INCDIR) -L$(LIBDIR) --f90flags="$(F95FLAGS)" \
		-c $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95\
		-l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS)
	mv _SHTOOLS.so pyshtools/.

pyshtools/_constant.so: $(SRCDIR)/PlanetsConstants.f95
	$(F2PY) --quiet --f90flags="$(F95FLAGS)" -c $(SRCDIR)/PlanetsConstants.f95 -m _constant
	mv _constant.so pyshtools/.

pyshtools/_SHTOOLS$(PY3EXT): $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95 fortran
	$(F2PY3) --quiet -I$(INCDIR) -L$(LIBDIR) --f90flags="$(F95FLAGS)" \
		-c $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95\
		-l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS)
	mv _SHTOOLS$(PY3EXT) pyshtools/.

pyshtools/_constant$(PY3EXT): $(SRCDIR)/PlanetsConstants.f95
	$(F2PY3) --quiet --f90flags="$(F95FLAGS)" -c $(SRCDIR)/PlanetsConstants.f95 -m _constant
	mv _constant$(PY3EXT) pyshtools/.

install: install-fortran

install-python2: python2
	mkdir -pv $(DESTDIR)$(SYSPYPATH)/pyshtools
	cp -R $(filter-out %$(PY3EXT), $(wildcard pyshtools/*)) $(DESTDIR)$(SYSPYPATH)/pyshtools/
	@echo
	@echo "Import shtools into Python 2 with:"
	@echo "----------------------------------"
	@echo "import sys"
	@echo "sys.path.append('$(SYSPYPATH)')"
	@echo "import pyshtools"
	@echo

install-python3: python3
	mkdir -pv $(DESTDIR)$(SYSPY3PATH)/pyshtools
	cp -R $(filter-out %.so, $(wildcard pyshtools/*)) $(DESTDIR)$(SYSPY3PATH)/pyshtools/
	cp -R $(filter %$(PY3EXT), $(wildcard pyshtools/*)) $(DESTDIR)$(SYSPY3PATH)/pyshtools/
	@echo
	@echo "Import shtools into Python 3 with:"
	@echo "----------------------------------"
	@echo "import sys"
	@echo "sys.path.append('$(SYSPY3PATH)')"
	@echo "import pyshtools"
	@echo

uninstall:
	-rm -r $(SYSPYPATH)/pyshtools/
	-rm -r $(SYSPY3PATH)/pyshtools/
	-rm -r $(SYSLIBPATH)/lib$(LIBNAME).a
	-rm -r $(SYSLIBPATH)/lib$(LIBNAMEMP).a
	-rm -r $(SYSMODPATH)/fftw3.mod
	-rm -r $(SYSMODPATH)/planetsconstants.mod
	-rm -r $(SYSMODPATH)/shtools.mod
	-rm -r $(SYSSHAREPATH)/shtools/examples/
	-rm -r $(SYSDOCPATH)/shtools/index.html
	-rm -r $(SYSDOCPATH)/shtools/www/
	$(MAKE) -C $(FDOCDIR) -f Makefile VERSION=$(VERSION) uninstall
	$(MAKE) -C $(PYDOCDIR) -f Makefile VERSION=$(VERSION) uninstall

install-fortran: fortran
	-rm -r $(DESTDIR)$(SYSMODPATH)/fftw3.mod
	-rm -r $(DESTDIR)$(SYSMODPATH)/planetsconstants.mod
	-rm -r $(DESTDIR)$(SYSMODPATH)/shtools.mod
	mkdir -pv $(DESTDIR)$(SYSLIBPATH)
	cp $(LIBDIR)/lib$(LIBNAME).a $(DESTDIR)$(SYSLIBPATH)/lib$(LIBNAME).a
	-cp $(LIBDIR)/lib$(LIBNAMEMP).a $(DESTDIR)$(SYSLIBPATH)/lib$(LIBNAMEMP).a
	mkdir -pv $(DESTDIR)$(SYSMODPATH)
	cp $(INCDIR)/*.mod $(DESTDIR)$(SYSMODPATH)/
	mkdir -pv $(DESTDIR)$(SYSSHAREPATH)/shtools
	cp -R examples $(DESTDIR)$(SYSSHAREPATH)/shtools/
	mkdir -pv $(DESTDIR)$(SYSSHAREPATH)/man/man1
	cp -R man/man1/ $(DESTDIR)$(SYSSHAREPATH)/man/man1/
	mkdir -pv $(DESTDIR)$(SYSDOCPATH)/shtools
	cp index.html $(DESTDIR)$(SYSDOCPATH)/shtools/index.html
	cp -R www $(DESTDIR)$(SYSDOCPATH)/shtools/
	awk '{gsub("../../lib","$(PREFIX)/lib");print}' "examples/fortran/Makefile" > "temp.txt"
	awk '{gsub("../../modules","$(PREFIX)/include");print}' "temp.txt" > "temp2.txt"
	cp temp2.txt "$(DESTDIR)$(SYSSHAREPATH)/shtools/examples/fortran/Makefile"
	rm temp.txt
	rm temp2.txt
	@echo
	@echo "Compile your Fortran code with the following flags:"
	@echo "---------------------------------------------------"
	@echo $(F95) $(SYSMODFLAG) $(F95FLAGS) -L$(SYSLIBPATH) -l$(LIBNAME) $(FFTW) -lm $(LAPACK) $(BLAS)
	@echo

doc:
	@$(MAKE) -C $(FDOCDIR) -f Makefile VERSION=$(VERSION)
	@$(MAKE) -C $(PYDOCDIR) -f Makefile VERSION=$(VERSION)
	@echo "--> Documentation created successfully"

remove-doc:
	@-rm -f man/man1/*.1
	@-rm -f doc/pages/mydoc/fdoc/*.md
	@-rm -f doc/pages/mydoc/pydoc/*.md
	@echo "--> Removed man files and web site source md files"

www:
	@cd $(WWWSRC) ; $(JEKYLL) build -d ../$(WWWDEST)

remove-www:
	@-rm -r $(WWWDEST)

notebooks2:
	@$(MAKE) -C $(NBDIR) -f Makefile JUPYTER=$(JUPYTER)
	@echo "--> Notebook html files created successfully with Python 2"

notebooks3:
	@$(MAKE) -C $(NBDIR) -f Makefile JUPYTER=$(JUPYTER3)
	@echo "--> Notebook html files created successfully with Python 3"

clean: clean-fortran-tests clean-python-tests clean-python2 clean-python3 clean-libs remove-www

clean-fortran-tests:
	@$(MAKE) -C $(FEXDIR) -f Makefile clean
	@echo "--> Removed Fortran test suite executables and files"

clean-python-tests:
	@$(MAKE) -C $(PEXDIR) -f Makefile clean
	@echo "--> Removed Python test suite executables and files"

clean-python2:
	@-rm -rf _SHTOOLS.so.dSYM/ _constant.so.dSYM/
	@-rm -rf pyshtools/_SHTOOLS.so.dSYM/ pyshtools/_constant.so.dSYM/
	@-rm -f *.so
	@-rm -f pyshtools/*.so
	@-rm -f pyshtools/*.pyc
	@-rm -rf pyshtools/__pycache__/
	@-rm -rf pyshtools/doc
	@echo "--> Removed Python 2 files"

clean-python3:
	@-rm -rf _SHTOOLS$(PY3EXT).dSYM/ _constant$(PY3EXT).dSYM/
	@-rm -rf pyshtools/_SHTOOLS$(PY3EXT).dSYM/ pyshtools/_constant$(PY3EXT).dSYM/
	@-rm -rf _SHTOOLS.so.dSYM/ _constant.so.dSYM/
	@-rm -rf pyshtools/_SHTOOLS.so.dSYM/ pyshtools/_constant.so.dSYM/
	@-rm -f *.so
	@-rm -f pyshtools/*.so
	@-rm -f pyshtools/*.pyc
	@-rm -rf pyshtools/__pycache__/
	@-rm -rf pyshtools/doc
	@echo "--> Removed Python 3 files"

clean-libs:
	@-$(MAKE) -C $(SRCDIR) -f Makefile clean
	@-rm -rf lib
	@-rm -rf modules
	@-rm -rf NONE
	@-rm -rf build
	@-rm -rf pyshtools.egg-info
	@-rm -f src/_SHTOOLS-f2pywrappers.f src/_SHTOOLSmodule.c
	@-rm -rf dist
	@echo "--> Removed lib, module, object files, compiled Python files and tests"
	@echo
	@echo \*\*\* If you installed pyshtools using \"pip install -e .\" you should
	@echo \*\*\* also execute \"pip uninstall pyshtools\".

remove-notebooks:
	@$(MAKE) -C $(NBDIR) -f Makefile clean
	@echo "--> Removed notebook html files"

fortran-tests: fortran
	@$(MAKE) -C $(FEXDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS)" LIBNAME="$(LIBNAME)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo "--> Make of Fortran test suite successful"
	@echo
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests
	@echo
	@echo "--> Ran all Fortran examples and tests"

fortran-tests-no-timing: fortran
	@$(MAKE) -C $(FEXDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS)" LIBNAME="$(LIBNAME)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo "--> Make of Fortran test suite successful"
	@echo
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests-no-timing
	@echo
	@echo "--> Ran all Fortran examples and tests"

fortran-tests-mp: fortran-mp
	@$(MAKE) -C $(FEXDIR) -f Makefile all F95=$(F95) F95FLAGS="$(OPENMPFLAGS) $(F95FLAGS)" LIBNAME="$(LIBNAMEMP)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo "--> Make of Fortran test suite successful"
	@echo
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests
	@echo
	@echo "--> Ran all Fortran examples and tests"

fortran-tests-no-timing-mp: fortran-mp
	@$(MAKE) -C $(FEXDIR) -f Makefile all F95=$(F95) F95FLAGS="$(OPENMPFLAGS) $(F95FLAGS)" LIBNAME="$(LIBNAMEMP)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo "--> Make of Fortran test suite successful"
	@echo
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests-no-timing
	@echo
	@echo "--> Ran all Fortran examples and tests"

run-fortran-tests: fortran
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests
	@echo
	@echo "--> Ran all Fortran examples and tests"

run-fortran-tests-no-timing: fortran
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests-no-timing
	@echo
	@echo "--> Ran all Fortran examples and tests"

python2-tests: python2
	@$(MAKE) -C $(PEXDIR) -f Makefile all PYTHON=$(PYTHON)
	@echo
	@echo "--> Ran all Python 2 tests"

python3-tests: python3
	@$(MAKE) -C $(PEXDIR) -f Makefile all PYTHON=$(PYTHON3)
	@echo
	@echo "--> Ran all Python 3 tests"

python2-tests-no-timing: python2
	@$(MAKE) -C $(PEXDIR) -f Makefile no-timing PYTHON=$(PYTHON)
	@echo
	@echo "--> Ran all Python 2 tests"

python3-tests-no-timing: python3
	@$(MAKE) -C $(PEXDIR) -f Makefile no-timing PYTHON=$(PYTHON3)
	@echo
	@echo "--> Ran all Python 3 tests"
