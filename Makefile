###############################################################################
#
#   INSTRUCTIONS
#
#   The normal user should only have to use the following commands
#
#       make (all, fortran)   : Install the Fortran 95 components.
#       make fortran-mp       : Install the Fortan OpenMP components.
#       make fortran-tests    : Compile and run the Fortran test/example suite.
#       make fortran-tests-mp : Compile and run the fortran test/example suite
#                               with OpenMp support.
#       make fortran-tests-no-timing    : Compile and run the Fortran
#                                         test/example suite (excluding timing
#                                         tests).
#       make fortran-tests-no-timing-mp : Compile and run the Fortran OpenMP
#                                         test/example suite (excluding timing
#                                         tests).
#       make python-tests     : Run the Python test/example suite.
#       make python-tests-no-timing  : Run the Python test/example suite,
#                                      (excluding timing tests).
#       make run-notebooks    : Run notebooks to test for errors.
#       make install          : Place the compiled libraries and docs in
#                               $(DESTDIR)$(PREFIX) [default is {}{usr/local}].
#       make uninstall        : Remove files copied to $(DESTDIR)$(PREFIX).
#       make clean            : Return the folder to its original state.
#       make check            : Check syntax of python files using flake8.
#
#
#   In some cases, where there are underscore problems when linking to the
#   LAPACK, BLAS, and FFTW3 libraries, it might be necessary to set:
#
#       LAPACK_UNDERSCORE=1
#
#
#   This Makefile accepts the following optional arguments that can be passed
#   using the syntax : make NAME="options"
#
#       F95         : Name and path of the Fortran-95 compiler.
#       F95FLAGS    : Fortran-95 compiler flags.
#       OPENMPFLAGS : Fortran-95 OpenMP compiler flags.
#       PYTHON      : Name (including path) of the Python executable.
#       FFTW        : Name and path of the FFTW3 library of the form
#                     "-Lpath -lfftw3 -lm".
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
#       Compile the Fortran 95 components.
#
#   make fortran-mp
#       Compile the fortran 95 component of the archive with OpenMP support.
#
#   make clean
#       Remove the compiled lib, module, object, and Python files. Also removes
#       compiled fortran and Python tests.
#
#   make run-fortran-tests
#       Run all Fortran examples and test suite.
#
#   make run-fortran-tests-no-timing
#       Run all Fortran examples and test suite (excluding the timing tests).
#
#   make run-fortran-tests-mp
#       Run all Fortran OpenMP examples and test suite.
#
#   make run-fortran-tests-no-timing-mp
#       Run all Fortran OpenMP examples and test suite (excluding the timing
#       tests).
#
#   make clean-fortran-tests
#       Delete compiled Fortran test suite programs.
#
#   make python-tests
#       Run all Python tests.
#
#   make python-tests-no-timing
#       Run all Python tests (with exception of timing/accuracy tests).
#
#   make clean-python-tests
#       Detele all compiled Python test files.
#
#   make clean-python
#       Detele all compiled Python files.
#
#   make notebooks
#       Convert notebooks html.
#
#   make remove-notebooks
#       Remove html notebooks.
#
#   make doc
#       Create the man pages from input markdown files and reformat
#       md files for use with the static web site. These files are PRE-MADE in
#       the distribution. To remake these files for a new release, it will be
#       necessary to install "pandoc", "ghc" and "cabal-install" (all using
#       brew on macOS), and then execute "cabal update" and
#       "cabal install --lib pandoc-types". If errors are encountered try
#       "brew uninstall pandoc" and "cabal install pandoc".
#
#   make remove-doc
#       Remove the man and html-man pages.
#
#   make www
#       Make the static html web documention in the directory www using Jekyll.
#       First, you must install "ruby" (using brew on macOS), and then install
#       the gem bundler using "gem install bundler jekyll". To serve the web
#       documents without making static html files, go to the directory `docs`
#       type the command `bundle update` and `bundle exec jekyll serve`, and
#       then open http://127.0.0.1:4000.
#
#   make remove-www
#       Remove the directory containing the static html web site.
#
###############################################################################

VERSION = 4.7

LIBNAME = SHTOOLS
LIBNAMEMP = SHTOOLS-mp

F95 = gfortran
PYTHON = python3
JUPYTER = jupyter nbconvert --ExecutePreprocessor.kernel_name=python3
JEKYLL = bundle exec jekyll
FLAKE8 = flake8
SHELL = /bin/sh
PY3EXT = $(shell $(PYTHON) -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))' || echo nopy3)


FDOCDIR = src/fdoc
PYDOCDIR = src/pydoc
SRCDIR = src
LIBDIR = lib
MODDIR = modules
FEXDIR = examples/fortran
PEXDIR = examples/python
NBDIR = examples/notebooks
WWWSRC = docs
WWWDEST = www
LIBPATH = $(PWD)/$(LIBDIR)
MODPATH = $(PWD)/$(MODDIR)
PYPATH = $(PWD)

PREFIX = /usr/local
SYSLIBPATH = $(PREFIX)/lib
SYSMODPATH = $(PREFIX)/include
SYSSHAREPATH = $(PREFIX)/share

FFTW = -L$(SYSLIBPATH) -lfftw3 -lm
LAPACK_UNDERSCORE = 0

FLAKE8_FILES = setup.py pyshtools examples/python

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
F95FLAGS ?= -m64 -fPIC -O3 -std=gnu -ffast-math
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


.PHONY: all fortran fortran-mp install fortran-tests fortran-tests-mp \
	fortran-tests-no-timing fortran-tests-no-timing-mp run-fortran-tests \
	run-fortran-tests-no-timing doc remove-doc python-tests \
	python-tests-no-timing uninstall clean clean-fortran-tests \
	clean-python-tests clean-python clean-libs remove-notebooks notebooks \
	run-notebooks www remove-www help check

help:
	@echo "Commands:"
	@echo ""
	@echo "  fortran                      Compile the Fortran 95 components"
	@echo "  fortran-mp                   Compile the Fortran 95 OpenMP components"
	@echo "  fortran-tests                Compile and run the Fortran 95 tests and examples"
	@echo "  fortran-tests-mp             Compile and run the Fortran 95 OpenMP tests and examples"
	@echo "  fortran-tests-no-timing      Do not run the timing tests"
	@echo "  fortran-tests-no-timing-mp   Do not run the timing tests"
	@echo "  python-tests                 Run the python tests and examples"
	@echo "  python-tests-no-timing       Do not run the timing tests"
	@echo "  run-notebooks                Execute the python notebooks"
	@echo "  clean                        Return the folder to its original state"
	@echo "  check                        Check syntax of python files using flake8"
	@echo ""

all: fortran

fortran:
	mkdir -pv $(LIBPATH)
	mkdir -pv $(MODPATH)
	$(MAKE) -C $(SRCDIR) -f Makefile all F95="$(F95)" F95FLAGS="$(F95FLAGS)" LIBNAME="$(LIBNAME)" LAPACK_FLAGS="$(LAPACK_FLAGS)" LIBPATH="$(LIBPATH)" MODPATH="$(MODPATH)"
	@echo "--> make fortran successful!"
	@echo
	@echo "Compile your Fortran code with the following flags:"
	@echo "---------------------------------------------------"
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAME) $(FFTW) $(LAPACK) $(BLAS)
	@echo

fortran-mp:
# Delete .o files before and after compiling with OpenMP to avoid issues with "fortran" build.
	mkdir -pv $(LIBPATH)
	mkdir -pv $(MODPATH)
	-$(MAKE) -C $(SRCDIR) -f Makefile clean-obs-mod
	$(MAKE) -C $(SRCDIR) -f Makefile all F95="$(F95)" F95FLAGS="$(OPENMPFLAGS) $(F95FLAGS)" LIBNAME="$(LIBNAMEMP)" LAPACK_FLAGS="$(LAPACK_FLAGS)" LIBPATH="$(LIBPATH)" MODPATH="$(MODPATH)"
	-$(MAKE) -C $(SRCDIR) -f Makefile clean-obs-mod
	@echo "--> make fortran-mp successful!"
	@echo
	@echo "Compile your Fortran code with the following flags:"
	@echo "---------------------------------------------------"
	@echo $(F95) $(MODFLAG) $(OPENMPFLAGS) $(F95FLAGS) -L$(LIBPATH) -l$(LIBNAMEMP) $(FFTW) $(LAPACK) $(BLAS)
	@echo

install: fortran
	mkdir -pv $(DESTDIR)$(SYSLIBPATH)
	cp $(LIBPATH)/lib$(LIBNAME).a $(DESTDIR)$(SYSLIBPATH)/lib$(LIBNAME).a
	-cp $(LIBPATH)/lib$(LIBNAMEMP).a $(DESTDIR)$(SYSLIBPATH)/lib$(LIBNAMEMP).a
	mkdir -pv $(DESTDIR)$(SYSMODPATH)
	cp $(MODPATH)/*.mod $(DESTDIR)$(SYSMODPATH)/
	mkdir -pv $(DESTDIR)$(SYSSHAREPATH)/examples/shtools
	cp -R examples/fortran/ $(DESTDIR)$(SYSSHAREPATH)/examples/shtools/
	mkdir -pv $(DESTDIR)$(SYSSHAREPATH)/man/man1
	cp -R man/man1/ $(DESTDIR)$(SYSSHAREPATH)/man/man1/
	awk '{gsub("../../lib","$(PREFIX)/lib");print}' "examples/fortran/Makefile" > "temp.txt"
	awk '{gsub("../../modules","$(PREFIX)/include");print}' "temp.txt" > "temp2.txt"
	cp temp2.txt "$(DESTDIR)$(SYSSHAREPATH)/examples/shtools/Makefile"
	rm temp.txt
	rm temp2.txt
	@echo
	@echo "Compile your Fortran code with the following flags:"
	@echo "---------------------------------------------------"
	@echo $(F95) $(SYSMODFLAG) $(F95FLAGS) -L$(SYSLIBPATH) -l$(LIBNAME) $(FFTW) $(LAPACK) $(BLAS)
	@echo

uninstall:
	-rm -r $(SYSLIBPATH)/lib$(LIBNAME).a
	-rm -r $(SYSLIBPATH)/lib$(LIBNAMEMP).a
	-rm -r $(SYSMODPATH)/fftw3.mod
	-rm -r $(SYSMODPATH)/planetsconstants.mod
	-rm -r $(SYSMODPATH)/shtools.mod
	-rm -r $(SYSMODPATH)/ftypes.mod
	-rm -r $(SYSSHAREPATH)/examples/shtools
	$(MAKE) -C $(FDOCDIR) -f Makefile uninstall PREFIX=$(PREFIX)

doc:
	@$(MAKE) -C $(FDOCDIR) -f Makefile VERSION=$(VERSION)
	@$(MAKE) -C $(PYDOCDIR) -f Makefile
	@echo "--> Documentation created successfully"

remove-doc:
	@-rm -f man/man1/*.1
	@-rm -f doc/pages/fortran/fdoc/*.md
	@-rm -f doc/pages/mydoc/pydoc/*.md
	@echo "--> Removed man files and web site source md files"

www:
	@cd $(WWWSRC) ; bundle update ; $(JEKYLL) build -d ../$(WWWDEST)

remove-www:
	@-rm -rf $(WWWDEST)

notebooks:
	@$(MAKE) -C $(NBDIR) -f Makefile JUPYTER="$(JUPYTER)"
	@echo "--> Notebook html files created successfully"

remove-notebooks:
	@$(MAKE) -C $(NBDIR) -f Makefile clean
	@echo "--> Removed notebook html files"

run-notebooks:
	@$(MAKE) -C $(NBDIR) -f Makefile run-notebooks
	@echo "--> Notebooks executed successfully"

clean: clean-fortran-tests clean-python-tests clean-python clean-libs remove-www

clean-fortran-tests:
	@$(MAKE) -C $(FEXDIR) -f Makefile clean
	@echo "--> Removed Fortran test suite executables and files"

clean-python-tests:
	@$(MAKE) -C $(PEXDIR) -f Makefile clean
	@echo "--> Removed Python test suite executables and files"

clean-python:
	@-rm -rf _SHTOOLS$(PY3EXT).dSYM/
	@-rm -rf pyshtools/_SHTOOLS$(PY3EXT).dSYM/
	@-rm -rf _SHTOOLS.so.dSYM/
	@-rm -rf pyshtools/_SHTOOLS.so.dSYM/
	@-rm -f *.so
	@-rm -rf __pycache__/
	@-rm -f pyshtools/*.so
	@-rm -f pyshtools/*.pyc
	@-rm -rf pyshtools/__pycache__/
	@-rm -rf pyshtools/doc
	@echo "--> Removed Python files"

clean-libs:
	@-$(MAKE) -C $(SRCDIR) -f Makefile clean
	@-rm -rf $(LIBPATH)
	@-rm -rf $(MODPATH)
	@-rm -rf NONE
	@-rm -rf build
	@-rm -rf pyshtools.egg-info
	@-rm -f src/_SHTOOLS-f2pywrappers.f src/_SHTOOLSmodule.c
	@-rm -rf dist
	@echo "--> Removed lib, module, object files, compiled Python files and tests"
	@echo
	@echo \*\*\* If you installed pyshtools using \"pip install -e .\" you should
	@echo \*\*\* also execute \"pip uninstall pyshtools\".

fortran-tests: fortran
	@$(MAKE) -C $(FEXDIR) -f Makefile all F95="$(F95)" F95FLAGS="$(F95FLAGS)" LIBNAME="$(LIBNAME)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo "--> Make of Fortran test suite successful"
	@echo
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests
	@echo
	@echo "--> Ran all Fortran examples and tests"

fortran-tests-no-timing: fortran
	@$(MAKE) -C $(FEXDIR) -f Makefile all F95="$(F95)" F95FLAGS="$(F95FLAGS)" LIBNAME="$(LIBNAME)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo "--> Make of Fortran test suite successful"
	@echo
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests-no-timing
	@echo
	@echo "--> Ran all Fortran examples and tests"

fortran-tests-mp: fortran-mp
	@$(MAKE) -C $(FEXDIR) -f Makefile all F95="$(F95)" F95FLAGS="$(OPENMPFLAGS) $(F95FLAGS)" LIBNAME="$(LIBNAMEMP)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo "--> Make of Fortran test suite successful"
	@echo
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests
	@echo
	@echo "--> Ran all Fortran examples and tests"

fortran-tests-no-timing-mp: fortran-mp
	@$(MAKE) -C $(FEXDIR) -f Makefile all F95="$(F95)" F95FLAGS="$(OPENMPFLAGS) $(F95FLAGS)" LIBNAME="$(LIBNAMEMP)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo "--> Make of Fortran test suite successful"
	@echo
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests-no-timing
	@echo
	@echo "--> Ran all Fortran examples and tests"

run-fortran-tests: fortran
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests F95="$(F95)" F95FLAGS="$(F95FLAGS)" LIBNAME="$(LIBNAME)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo "--> Ran all Fortran examples and tests"

run-fortran-tests-no-timing: fortran
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests-no-timing F95="$(F95)" F95FLAGS="$(F95FLAGS)" LIBNAME="$(LIBNAME)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo "--> Ran all Fortran examples and tests"

run-fortran-tests-mp: fortran-mp
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests LIBNAME="$(LIBNAMEMP)" F95="$(F95)" F95FLAGS="$(OPENMPFLAGS) $(F95FLAGS)" LIBNAME="$(LIBNAMEMP)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo "--> Ran all Fortran examples and tests"

run-fortran-tests-no-timing-mp: fortran-mp
	@$(MAKE) -C $(FEXDIR) -f Makefile run-fortran-tests-no-timing LIBNAME="$(LIBNAMEMP)" F95="$(F95)" F95FLAGS="$(OPENMPFLAGS) $(F95FLAGS)" LIBNAME="$(LIBNAMEMP)" FFTW="$(FFTW)" LAPACK="$(LAPACK)" BLAS="$(BLAS)"
	@echo
	@echo "--> Ran all Fortran examples and tests"

python-tests:
	@$(MAKE) -C $(PEXDIR) -f Makefile all PYTHON=$(PYTHON)
	@echo
	@echo "--> Ran all Python tests"

python-tests-no-timing:
	@$(MAKE) -C $(PEXDIR) -f Makefile no-timing PYTHON=$(PYTHON)
	@echo
	@echo "--> Ran all Python tests"

check:
	@$(FLAKE8) --extend-ignore=E741,W605 --exclude=versioneer.py,pyshtools/_version.py $(FLAKE8_FILES)
