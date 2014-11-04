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
#	make python-wrapper
#		Compile the python wrapper with the f2py compiler. This should be 
#		after the Fortran files are compiled.
#
#	make clean
#		Remove the compile lib, module, and object files.
#
#	make examples
#		Compile example programs. Optionally, one
#		can specify the parameters F95="my compiler" and 
#		F95FLAGS="my compiler flags", which should be identical to
#		those used to make "all".
#
#	make remove-examples
#		Delete compiled example programs.
#
#	make run-examples
#		Run all examples programs.
#
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

VERSION = 2.10

F95 = gfortran
F2PY = f2py

SHELL  = /bin/tcsh
MAKE   = make
DOCDIR = src/doc
SRCDIR = src
LIBDIR = lib
INCDIR = modules
EXDIR  = examples


.PHONY: all all2 all3 python-wrapper install doc remove-doc clean getflags examples remove-examples run-examples
	
all: getflags
	$(MAKE) -C $(SRCDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE SUCCESSFUL!
	@echo	
	@echo ---------------------------------------------------------------------------------------------------
	@echo Compile your code with the following flags:
	@echo
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -Llibpath -lSHTOOLS -lfftw3 -lm
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
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -Llibpath -lSHTOOLS -lfftw3 -lm
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
	@echo $(F95) $(MODFLAG) $(F95FLAGS) -Llibpath -lSHTOOLS -lfftw3 -lm
	@echo
	@echo where modpath and libpath are replaced with their respective paths.
	@echo ---------------------------------------------------------------------------------------------------
	@echo 

python-wrapper: all
	$(F2PY) -I$(INCDIR) -L$(LIBDIR) --f90flags="-m64 -fPIC" \
	    -c $(SRCDIR)/pyshtools.pyf $(SRCDIR)/PythonWrapper.f95\
	    -lSHTOOLS -lfftw3 -lm -llapack -lblas	
	$(F2PY) -c $(SRCDIR)/PlanetsConstants.f95 -m constants
	mv _SHTOOLS.so pyshtools/.
	@echo
	@echo MAKE SUCCESSFUL!
	@echo	
	@echo ---------------------------------------------------------------------------------------------------
	@echo import shtools into Python with:
	@echo
	@echo import shtools
	@echo ---------------------------------------------------------------------------------------------------
	@echo 
	
getflags:
ifeq ($(F95),f95)
# Default Absoft Pro Fortran flags
F95FLAGS ?= -m64 -O3 -YEXT_NAMES=LCS -YEXT_SFX=_ -fpic #-march=host
MODFLAG = -p modpath
endif

ifeq ($(F95),gfortran)
# Default gfortran flags
F95FLAGS ?= -m64 -fPIC -O3 # -march=native
MODFLAG = -Imodpath
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
	$(MAKE) -C $(DOCDIR) -f Makefile VERSION=$(VERSION)
	@echo DOCUMENTATION SUCCESSFULLY CREATED

remove-doc:
	@rm -f man/man1/*.1
	@rm -f www/man/*.html
	@echo
	@echo REMOVED MAN AND HTML-MAN FILES
	
clean:
	$(MAKE) -C $(SRCDIR) -f Makefile clean
	@echo
	@echo REMOVED LIB, MODULE, AND OBJECT FILES

examples: getflags
	$(MAKE) -C $(EXDIR) -f Makefile all F95=$(F95) F95FLAGS="$(F95FLAGS)"
	@echo
	@echo MAKE OF EXAMPLES SUCCESSFUL

remove-examples:
	$(MAKE) -C $(EXDIR) -f Makefile clean
	@echo
	@echo REMOVED EXAMPLE EXECUTABLES

run-examples:
	$(MAKE) -C $(EXDIR) -f Makefile run-examples
	@echo
	@echo RAN ALL EXAMPLE PROGRAMS

	
