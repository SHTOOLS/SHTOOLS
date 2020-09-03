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
	SHBiasAdmitCorr.f95 SHCimToVector.f95 YilmIndexVector.f95 \
	ComputeDMap.f95 SHReturnTapersMap.f95 Curve2Mask.f95 MakeEllipseCoord.f95 \
	CilmPlusRhoH.f95 CilmMinusRhoH.f95 BAtoHilmRhoH.f95 NormalGravity.f95 \
	SHMultiTaperMaskSE.f95 SHMultiTaperMaskCSE.f95 SHBiasKMask.f95 \
	SHRotateTapers.f95 SlepianCoeffs.f95 SlepianCoeffsToSH.f95 \
	SHSCouplingMatrix.f95 SHMTVar.f95 SHSlepianVar.f95 SHSCouplingMatrixCap.f95 \
	cWrapper.f95

OBJS0 = CilmPlus.o CilmMinus.o ComputeDG82.o \
	ComputeDm.o DHaj.o djpi2.o BAtoHilm.o MakeGrid2D.o GLQGridCoord.o \
	MakeGridPoint.o MakeGridPointC.o PlBar.o PlBar_d1.o \
	PLegendre.o PLegendre_d1.o PLegendreA.o PLegendreA_d1.o PlmBar.o \
	PlmBar_d1.o PlmIndex.o PlmSchmidt.o PlmSchmidt_d1.o PlSchmidt.o \
	PlSchmidt_d1.o PreGLQ.o Random.o SHAdmitCorr.o SHBias.o SHBiasK.o \
	SHConvertCoef.o SHFindLWin.o SHGLQ.o SHLocalizedAdmitCorr.o \
	SHMultiply.o SHMultiTaperCSE.o SHMultiTaperSE.o SHPowerSpectra.o \
	SHRead.o SHReadJPL.o SHReturnTapersM.o SHRotateCoef.o \
	SHRotateRealCoef.o SphericalCapCoef.o Wigner3j.o \
	DownContFilter.o SHRead2.o MakeGeoidGrid.o \
	MakeCircleCoord.o SHMTCouplingMatrix.o SHReturnTapers.o SHSjkPG.o \
	PlON.o PlON_d1.o PlmON.o PlmON_d1.o \
	SHConfidence.o SHMagPowerSpectra.o \
	SHPowerSpectraC.o SHBiasAdmitCorr.o SHCilmToVector.o \
	YilmIndexVector.o ComputeDMap.o SHReturnTapersMap.o Curve2Mask.o \
	MakeEllipseCoord.o CilmPlusRhoH.o CilmMinusRhoH.o BAtoHilmRhoH.o \
	NormalGravity.o SHMultiTaperMaskSE.o SHMultiTaperMaskCSE.o SHBiasKMask.o \
	SHRotateTapers.o SlepianCoeffs.o SlepianCoeffsToSH.o SHSCouplingMatrix.o \
	SHMTVar.o SHSlepianVar.o SHSCouplingMatrixCap.o \
	cWrapper.o

SRCSLAPACK = EigValSym.F95 EigValVecSym.F95 EigValVecSymTri.F95 \
	SHExpandLSQ.F95 SHMTDebias.F95 SHMTVarOpt.F95

OBJSLAPACK = EigValSym.o EigValVecSym.o EigValVecSymTri.o SHExpandLSQ.o \
	SHMTDebias.o SHMTVarOpt.o

SRCSFFTW = MakeGridDH.f95 MakeGridDHC.f95 MakeGridGLQ.f95 MakeGridGLQC.f95 \
	SHExpandDH.f95 SHExpandDHC.f95 SHExpandGLQ.f95 SHExpandGLQC.f95 \
	MakeGravGradGridDH.f95 MakeGravGridDH.f95 MakeMagGridDH.f95 \
	MakeMagGradGridDH.f95

OBJSFFTW = MakeGridDH.o MakeGridDHC.o MakeGridGLQ.o MakeGridGLQC.o \
	SHExpandDH.o SHExpandDHC.o SHExpandGLQ.o SHExpandGLQC.o \
	MakeGravGradGridDH.o MakeGravGridDH.o MakeMagGridDH.o \
	MakeMagGradGridDH.o

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
%.o: %.f95
	$(F95) -c $(F95FLAGS) $<
%.o: %.F95
	$(F95) -c $(F95FLAGS) $(LAPACK_FLAGS) $<


ftypes.mod : ftypes.o
SHTOOLS.mod : ftypes.mod
FFTW3.mod : ftypes.mod
PlanetsConstants.mod : ftypes.mod

cWrapper.o : SHTOOLS.mod
BAtoHilm.o : ftypes.mod SHTOOLS.mod
BAtoHilmRhoH.o : ftypes.mod SHTOOLS.mod
CilmMinus.o : ftypes.mod SHTOOLS.mod
CilmMinusRhoH.o : ftypes.mod SHTOOLS.mod
CilmPlus.o : ftypes.mod SHTOOLS.mod
CilmPlusRhoH.o : ftypes.mod SHTOOLS.mod
ComputeDG82.o : ftypes.mod
ComputeDMap.o : ftypes.mod SHTOOLS.mod
ComputeDm.o : ftypes.mod SHTOOLS.mod
Curve2Mask.o : ftypes.mod
DHaj.o : ftypes.mod
DownContFilter.o : ftypes.mod
EigValSym.o : ftypes.mod
EigValVecSym.o : ftypes.mod
EigValVecSymTri.o : ftypes.mod
GLQGridCoord.o : ftypes.mod SHTOOLS.mod
MakeCircleCoord.o : ftypes.mod
MakeEllipseCoord.o : ftypes.mod
MakeGeoidGrid.o : ftypes.mod SHTOOLS.mod
MakeGravGradGridDH.o : ftypes.mod FFTW3.mod
MakeGravGridDH.o : ftypes.mod FFTW3.mod SHTOOLS.mod
MakeGrid2D.o : ftypes.mod SHTOOLS.mod
MakeGridDH.o : ftypes.mod FFTW3.mod SHTOOLS.mod
MakeGridDHC.o : ftypes.mod FFTW3.mod SHTOOLS.mod
MakeGridGLQ.o : ftypes.mod FFTW3.mod SHTOOLS.mod
MakeGridGLQC.o : ftypes.mod FFTW3.mod SHTOOLS.mod
MakeGridPoint.o : ftypes.mod SHTOOLS.mod
MakeGridPointC.o : ftypes.mod SHTOOLS.mod
MakeMagGradGridDH.o : ftypes.mod FFTW3.mod
MakeMagGridDH.o : ftypes.mod FFTW3.mod SHTOOLS.mod
NormalGravity.o : ftypes.mod
PLegendre.o : ftypes.mod
PLegendreA.o : ftypes.mod SHTOOLS.mod
PLegendreA_d1.o : ftypes.mod SHTOOLS.mod
PLegendre_d1.o : ftypes.mod
PlBar.o : ftypes.mod
PlBar_d1.o : ftypes.mod
PlON.o : ftypes.mod
PlON_d1.o : ftypes.mod
PlSchmidt.o : ftypes.mod
PlSchmidt_d1.o : ftypes.mod
PlmBar.o : ftypes.mod SHTOOLS.mod
PlmBar_d1.o : ftypes.mod SHTOOLS.mod
PlmON.o : ftypes.mod SHTOOLS.mod
PlmON_d1.o : ftypes.mod SHTOOLS.mod
PlmSchmidt.o : ftypes.mod SHTOOLS.mod
PlmSchmidt_d1.o : ftypes.mod SHTOOLS.mod
PreGLQ.o : ftypes.mod
Random.o : ftypes.mod SHTOOLS.mod
SHAdmitCorr.o : ftypes.mod SHTOOLS.mod
SHBias.o : ftypes.mod SHTOOLS.mod
SHBiasAdmitCorr.o : ftypes.mod SHTOOLS.mod
SHBiasK.o : ftypes.mod SHTOOLS.mod
SHBiasKMask.o : ftypes.mod SHTOOLS.mod
SHCilmToVector.o : ftypes.mod
SHConfidence.o : ftypes.mod
SHConvertCoef.o : ftypes.mod
SHExpandDH.o : ftypes.mod FFTW3.mod SHTOOLS.mod
SHExpandDHC.o : ftypes.mod FFTW3.mod SHTOOLS.mod
SHExpandGLQ.o : ftypes.mod FFTW3.mod SHTOOLS.mod
SHExpandGLQC.o : ftypes.mod FFTW3.mod SHTOOLS.mod
SHExpandLSQ.o : ftypes.mod SHTOOLS.mod
SHFindLWin.o : ftypes.mod SHTOOLS.mod
SHGLQ.o : ftypes.mod SHTOOLS.mod
SHLocalizedAdmitCorr.o : ftypes.mod SHTOOLS.mod
SHMTCouplingMatrix.o : ftypes.mod SHTOOLS.mod
SHMTDebias.o : ftypes.mod SHTOOLS.mod
SHMTVar.o : ftypes.mod  SHTOOLS.mod
SHMTVarOpt.o : ftypes.mod SHTOOLS.mod
SHMagPowerSpectra.o : ftypes.mod
SHMultiTaperCSE.o : ftypes.mod SHTOOLS.mod
SHMultiTaperMaskCSE.o : ftypes.mod SHTOOLS.mod
SHMultiTaperMaskSE.o : ftypes.mod SHTOOLS.mod
SHMultiTaperSE.o : ftypes.mod SHTOOLS.mod
SHMultiply.o : ftypes.mod SHTOOLS.mod
SHPowerSpectra.o : ftypes.mod
SHPowerSpectraC.o : ftypes.mod
SHRead.o : ftypes.mod
SHRead2.o : ftypes.mod
SHReadJPL.o : ftypes.mod
SHReturnTapers.o : ftypes.mod SHTOOLS.mod
SHReturnTapersM.o : ftypes.mod SHTOOLS.mod
SHReturnTapersMap.o : SHTOOLS.mod
SHRotateCoef.o : ftypes.mod
SHRotateRealCoef.o : ftypes.mod SHTOOLS.mod
SHRotateTapers.o : ftypes.mod SHTOOLS.mod
SHSCouplingMatrix.o : ftypes.mod
SHSCouplingMatrixCap.o : ftypes.mod
SHSjkPG.o : ftypes.mod SHTOOLS.mod
SHSlepianVar.o : ftypes.mod
SlepianCoeffs.o: ftypes.mod SHTOOLS.mod
SlepianCoeffsToSH.o : ftypes.mod SHTOOLS.mod
SphericalCapCoef.o : ftypes.mod SHTOOLS.mod
Wigner3j.o : ftypes.mod
djpi2.o : ftypes.mod
