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

SHELL=/bin/sh

LIBTOOL = libtool
LIBTOOLFLAGS = -static
AR = ar
ARFLAGS = -r
RLIB = ranlib
RLIBFLAGS =
LDFLAGS = -s

LIBNAME = SHTOOLS
PROG = ../lib/lib$(LIBNAME).a

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

SRCS0 = ftypes.f95 SHTOOLS.f95 FFTW3.f95 CilmPlus.f95 CilmMinus.f95 \
	ComputeDG82.f95 ComputeDm.f95 DHaj.f95 djpi2.f95 BAtoHilm.f95 \
	MakeGrid2D.f95 GLQGridCoord.f95 MakeGridPoint.f95 MakeGridPointC.f95 \
	PlanetsConstants.f95 PlBar.f95 PlBar_d1.f95 PLegendre.f95 \
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
	SHSCouplingMatrix.f95 SHMTVar.f95 SHSlepianVar.f95 SHSCouplingMatrixCap.f95

OBJS0 = ftypes.o SHTOOLS.o FFTW3.o CilmPlus.o CilmMinus.o ComputeDG82.o \
	ComputeDm.o DHaj.o djpi2.o BAtoHilm.o MakeGrid2D.o GLQGridCoord.o \
	MakeGridPoint.o MakeGridPointC.o PlanetsConstants.o PlBar.o PlBar_d1.o \
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
	SHMTVar.o SHSlepianVar.o SHSCouplingMatrixCap.o

SRCSLAPACK = EigValSym.F95 EigValVecSym.F95 EigValVecSymTri.F95 \
	SHExpandLSQ.F95 SHMTDebias.F95 SHMTVarOpt.F95

OBJSLAPACK = EigValSym.o EigValVecSym.o EigValVecSymTri.o SHExpandLSQ.o \
	SHMTDebias.o SHMTVarOpt.o

SRCSFFTW = MakeGridDH.F95 MakeGridDHC.F95 MakeGridGLQ.F95 MakeGridGLQC.F95 \
	SHExpandDH.F95 SHExpandDHC.F95 SHExpandGLQ.F95 SHExpandGLQC.F95 \
	MakeGravGradGridDH.F95 MakeGravGridDH.F95 MakeMagGridDH.F95 \
	MakeMagGradGridDH.F95

OBJSFFTW = MakeGridDH.o MakeGridDHC.o MakeGridGLQ.o MakeGridGLQC.o \
	SHExpandDH.o SHExpandDHC.o SHExpandGLQ.o SHExpandGLQC.o \
	MakeGravGradGridDH.o MakeGravGridDH.o MakeMagGridDH.o \
	MakeMagGradGridDH.o

SRCS = $(SRCS0) $(SRCSLAPACK) $(SRCSFFTW)
OBJS = $(OBJS0) $(OBJSLAPACK) $(OBJSFFTW)


all: $(PROG)

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
#	@mv -f *.mod ../modules/
	@cp -f *.mod ../modules/
	@echo "--> Library and module files moved to ../lib and ../modules"
	@echo "--> Archive created successfully"

.PHONY: clean clean-obs-mod clean-prog-mod all

clean: clean-obs-mod clean-prog-mod

clean-obs-mod:
	@-rm -f $(OBJS0)
	@-rm -f $(OBJSLAPACK)
	@-rm -f $(OBJSFFTW)
	@-rm -f *.mod

clean-prog-mod:
	@-rm -f $(PROG)
	@-rm -f ../modules/*.mod
	
%.o: %.f95
	$(F95) -c $(F95FLAGS) $<
%.o: %.F95
	$(F95) -c $(F95FLAGS) $(LAPACK_FLAGS) $(FFTW3_FLAGS) $<


BAtoHilm.o: ftypes.o SHTOOLS.o
BAtoHilmRhoH.o : ftypes.o SHTOOLS.o
CilmMinus.o: ftypes.o SHTOOLS.o
CilmMinusRhoH.o : ftypes.o SHTOOLS.o
CilmPlus.o: ftypes.o SHTOOLS.o
CilmPlusRhoH.o : ftypes.o SHTOOLS.o
ComputeDG82.o : ftypes.o
ComputeDMap.o : ftypes.o SHTOOLS.o
ComputeDm.o: ftypes.o SHTOOLS.o
Curve2Mask.o : ftypes.o
DHaj.o : ftypes.o
DownContFilter.o : ftypes.o
EigValSym.o : ftypes.o
EigValVecSym.o : ftypes.o
EigValVecSymTri.o : ftypes.o
GLQGridCoord.o : ftypes.o SHTOOLS.o
MakeCircleCoord.o : ftypes.o
MakeEllipseCoord.o : ftypes.o
MakeGeoidGrid.o : ftypes.o SHTOOLS.o
MakeGravGradGridDH.o : ftypes.o FFTW3.o
MakeGravGridDH.o : ftypes.o FFTW3.o SHTOOLS.o
MakeGrid2D.o : ftypes.o SHTOOLS.o
MakeGridDH.o : ftypes.o FFTW3.o SHTOOLS.o
MakeGridDHC.o : ftypes.o FFTW3.o SHTOOLS.o
MakeGridGLQ.o : ftypes.o FFTW3.o SHTOOLS.o
MakeGridGLQC.o : ftypes.o FFTW3.o SHTOOLS.o
MakeGridPoint.o : ftypes.o SHTOOLS.o
MakeGridPointC.o : ftypes.o SHTOOLS.o
MakeMagGradGridDH.o : ftypes.o FFTW3.o
MakeMagGridDH.o : ftypes.o FFTW3.o SHTOOLS.o
NormalGravity.o : ftypes.o
PLegendre.o : ftypes.o
PLegendreA.o : ftypes.o SHTOOLS.o
PLegendreA_d1.o : ftypes.o SHTOOLS.o
Plegendre_d1.o : ftypes.o
PlBar.o : ftypes.o
PlBar_d1.o : ftypes.o
PlON.o : ftypes.o
PlON_d1.o : ftypes.o
PlSchmidt.o : ftypes.o
PlSchmidt_d1.o : ftypes.o
PlanetsConstants.o : ftypes.o
PlmBar.o : ftypes.o SHTOOLS.o
PlmBar_d1.o : ftypes.o SHTOOLS.o
PlmON.o : ftypes.o SHTOOLS.o
PlmON_d1.o : ftypes.o SHTOOLS.o
PlmSchmidt.o : ftypes.o SHTOOLS.o
PlmSchmidt_d1.o : ftypes.o SHTOOLS.o
PreGLQ.o : ftypes.o
Random.o : ftypes.o
SHAdmitCorr.o : ftypes.o SHTOOLS.o
SHBias.o : ftypes.o SHTOOLS.o
SHBiasAdmitCorr.o : ftypes.o SHTOOLS.o
SHBiasK.o : ftypes.o SHTOOLS.o
SHBiasKMask.o : ftypes.o SHTOOLS.o
SHCilmToVector.o : ftypes.o
SHConfidence.o : ftypes.o
SHConvertCoef.o : ftypes.o
SHExpandDH.o : ftypes.o FFTW3.o SHTOOLS.o
SHExpandDHC.o : ftypes.o FFTW3.o SHTOOLS.o
SHExpandGLQ.o : ftypes.o FFTW3.o SHTOOLS.o
SHExpandGLQC.o : ftypes.o FFTW3.o SHTOOLS.o
SHExpandLSQ.o : ftypes.o SHTOOLS.o
SHFindLWin.o : ftypes.o SHTOOLS.o
SHGLQ.o : ftypes.o SHTOOLS.o
SHLocalizedAdmitCorr.o : ftypes.o SHTOOLS.o
SHMTCouplingMatrix.o : ftypes.o SHTOOLS.o
SHMTDebias.o : ftypes.o SHTOOLS.o
SHMTVar.o : ftypes.o  SHTOOLS.o
SHMTVarOpt.o : ftypes.o SHTOOLS.o
SHMagPowerSpectra.o : ftypes.o
SHMultiTaperCSE.o : ftypes.o SHTOOLS.o
SHMultiTaperMaskCSE.o : ftypes.o SHTOOLS.o
SHMultiTaperMaskSE.o : ftypes.o SHTOOLS.o
SHMultiTaperSE.o : ftypes.o SHTOOLS.o
SHMultiply.o : ftypes.o SHTOOLS.o
SHPowerSpectra.o : ftypes.o
SHPowerSpectraC.o : ftypes.o
SHRead.o : ftypes.o
SHRead2.o : ftypes.o
SHReadJPL.o : ftypes.o
SHReturnTapers.o : ftypes.o SHTOOLS.o
SHReturnTapersM.o : ftypes.o SHTOOLS.o
SHReturnTapersMap.o : SHTOOLS.o
SHRotateCoef.o : ftypes.o
SHRotateRealCoef.o : ftypes.o SHTOOLS.o
SHRotateTapers.o : ftypes.o SHTOOLS.o
SHSCouplingMatrix.o : ftypes.o
SHSCouplingMatrixCap.o : ftypes.o
SHSjkPG.o : ftypes.o SHTOOLS.o
SHSlepianVar.o : ftypes.o
SlepianCoeffs.o: ftypes.o SHTOOLS.o
SlepianCoeffsToSH.o : ftypes.o SHTOOLS.o
SphericalCapCoef.o : ftypes.o SHTOOLS.o
Wigner3j.o : ftypes.o
djpi2.o : ftypes.o
