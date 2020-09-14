###############################################################################
#
#   make all
#       Compile program in the current directory.
#
#   make clean
#       Remove the lib and module files in the LOCAL directory.
#
###############################################################################

F95 = gfortran

SHELL = /bin/sh

LIBNAME = SHTOOLS
LIBPATH = ../lib
MODPATH = ../modules
PROG = $(LIBPATH)/lib$(LIBNAME).a

LIBTOOL = libtool
LIBTOOLFLAGS = -static
AR = ar
ARFLAGS = -r
RLIB = ranlib
RLIBFLAGS =
LDFLAGS = -s

ifeq ($(F95),f95)
# Default Absoft Pro Fortran flags
F95FLAGS ?= -m64 -O3 -YEXT_NAMES=LCS -YEXT_SFX=_ -fpic -speed_math=11 # -march=host
else ifeq ($(F95),g95)
# Default g95 flags.
F95FLAGS ?= -O3 -fno-second-underscore
else ifeq ($(F95),gfortran)
# Default gfortran flags
F95FLAGS ?= -m64 -fPIC -O3 -ffast-math
else ifeq ($(F95),ifort)
# Default intel fortran flags
F95FLAGS ?= -m64 -fpp -free -O3 -Tf
else ifeq ($(origin F95FLAGS), undefined)
F95FLAGS = -m64 -O3
endif

SRCSMOD = ftypes.f95 SHTOOLS.f95 FFTW3.f95 PlanetsConstants.f95

MOD = ftypes.mod SHTOOLS.mod FFTW3.mod PlanetsConstants.mod

OBJSMOD = ftypes.o SHTOOLS.o FFTW3.o PlanetsConstants.o

SRCS0 = CilmPlus.f95 CilmMinus.f95 \
	ComputeDG82.f95 ComputeDm.f95 DHaj.f95 djpi2.f95 BAtoHilm.f95 \
	MakeGrid2D.f95 GLQGridCoord.f95 MakeGridPoint.f95 MakeGridPointC.f95 \
	PlBar.f95 PlBar_d1.f95 PLegendre.f95 \
	PLegendre_d1.f95 PLegendreA.f95 PLegendreA_d1.f95 PlmBar.f95 \
	PlmBar_d1.f95 PlmIndex.f95 PlmSchmidt.f95 PlmSchmidt_d1.f95 \
	PlSchmidt.f95 PlSchmidt_d1.f95 PreGLQ.f95 Random.f95 \
	SHAdmitCorr.f95 SHBias.f95 SHBiasK.f95 SHConvertCoef.f95 \
	SHFindLWin.f95 SHGLQ.f95 SHLocalizedAdmitCorr.f95 SHMultiply.f95 \
	SHMultiTaperCSE.f95 SHMultiTaperSE.f95 SHPowerSpectra.f95 SHRead.f95 \
	SHReadJPL.f95 SHReturnTapersM.f95 SHRotateCoef.f95 \
	SHRotateRealCoef.f95 SphericalCapCoef.f95 \
	Wigner3j.f95 DownContFilter.f95 SHRead2.f95 MakeGeoidGrid.f95 \
	MakeCircleCoord.f95 SHMTCouplingMatrix.f95 \
	SHReturnTapers.f95 SHSjkPG.f95 PlON.f95 PlON_d1.f95 PlmON.f95 \
	PlmON_d1.f95 SHConfidence.f95 SHMagPowerSpectra.f95 SHPowerSpectraC.f95 \
	SHBiasAdmitCorr.f95 SHCilmToVector.f95 YilmIndexVector.f95 \
	ComputeDMap.f95 SHReturnTapersMap.f95 Curve2Mask.f95 MakeEllipseCoord.f95 \
	CilmPlusRhoH.f95 CilmMinusRhoH.f95 BAtoHilmRhoH.f95 NormalGravity.f95 \
	SHMultiTaperMaskSE.f95 SHMultiTaperMaskCSE.f95 SHBiasKMask.f95 \
	SHRotateTapers.f95 SlepianCoeffs.f95 SlepianCoeffsToSH.f95 \
	SHSCouplingMatrix.f95 SHMTVar.f95 SHSlepianVar.f95 SHSCouplingMatrixCap.f95

OBJS0 := $(SRCS0:%.f95=%.o)


SRCSLAPACK = EigValSym.F95 EigValVecSym.F95 EigValVecSymTri.F95 \
	SHExpandLSQ.F95 SHMTDebias.F95 SHMTVarOpt.F95

OBJSLAPACK := $(SRCSLAPACK:%.F95=%.o)

SRCSFFTW = MakeGridDH.f95 MakeGridDHC.f95 MakeGridGLQ.f95 MakeGridGLQC.f95 \
	SHExpandDH.f95 SHExpandDHC.f95 SHExpandGLQ.f95 SHExpandGLQC.f95 \
	MakeGravGradGridDH.f95 MakeGravGridDH.f95 MakeMagGridDH.f95 \
	MakeMagGradGridDH.f95

OBJSFFTW := $(SRCSFFTW:%.f95=%.o)

SRCS = $(SRCSMOD) $(SRCS0) $(SRCSLAPACK) $(SRCSFFTW)
OBJS = $(OBJSMOD) $(OBJS0) $(OBJSLAPACK) $(OBJSFFTW)

all: $(MOD) $(PROG)

$(PROG): $(OBJS)
	@echo
	@echo "--> Compilation of source files successful"
	@echo
	@rm -f $(PROG)
#	@$(LIBTOOL) $(LIBTOOLFLAGS) -o $(PROG) $(OBJS)
#	If you prefer to use libtool, uncomment the above line, and comment the two lines below (AR and RLIB)
	$(AR) $(ARFLAGS) $(PROG) $(OBJS)
	$(RLIB) $(RLIBFLAGS) $(PROG)
	@echo
	@echo "--> Creation of static library successful"
#	@rm -f $(OBJS)
#	@mv -f *.mod $(MODPATH)/
	@cp -f *.mod $(MODPATH)/
	@echo "--> Library and module files moved to $(LIBPATH) and $(MODPATH)"
	@echo "--> Archive created successfully"

.PHONY: clean clean-obs-mod clean-prog-mod all

clean: clean-obs-mod clean-prog-mod

clean-obs-mod:
	@-rm -f $(OBJSMOD)
	@-rm -f $(OBJS0)
	@-rm -f $(OBJSMOD)
	@-rm -f $(OBJSLAPACK)
	@-rm -f $(OBJSFFTW)
	@-rm -f *.mod

clean-prog-mod:
	@-rm -f $(PROG)
	@-rm -f $(MODPATH)/*.mod
	
%.mod: %.f95
	$(F95) -c $(F95FLAGS) $<
%.o: %.f95 SHTOOLS.mod
	$(F95) -c $(F95FLAGS) $<
%.o: %.F95 SHTOOLS.mod
	$(F95) -c $(F95FLAGS) $(LAPACK_FLAGS) $<

ftypes.mod : ftypes.f95
SHTOOLS.mod : ftypes.mod
FFTW3.mod : ftypes.mod
PlanetsConstants.mod : ftypes.mod
