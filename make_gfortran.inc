###############################################################################
##########  GNU Fortran Compiler  #############################################

FCHOME = /usr
FC = $(FCHOME)/bin/gfortran
COMPINC = -I./ -I$(FCHOME)/include -I$(FCHOME)/lib64/gfortran/modules -I$(INSTALLDIR)/include
COMPLIBS = -L./ -L$(FCHOME)/lib64 -L${INSTALLDIR}/lib

LIBS = $(COMPLIBS) $(COMPINC)

# Debugging flags
ifeq ($(RUN), DEBUG)
    # consider including -lefence for debugging
    FFLAGS = -O0 -g3 -Wall -Wextra -fimplicit-none -Wline-truncation -Wcharacter-truncation -Waliasing -Wimplicit-interface -Wunused-parameter -fwhole-file -fcheck=all -std=f2008 -pedantic -fbacktrace -Wunderflow -ffpe-trap=invalid,zero,overflow -Wuninitialized -Wmaybe-uninitialized -ffpe-summary=all -finit-integer=-2147483647 -finit-real=snan -fno-unsafe-math-optimizations -fsignaling-nans -frounding-math -Wno-surprising
endif
# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g -pg -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math
endif
# Production run flags
ifeq ($(RUN), OPT)
    FFLAGS = -O3 -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -ftree-vectorize -funroll-loops -finline-functions -flto=3
endif
  # Preprocessing flags
FPPFLAGS = -x f95-cpp-input $(ncFPPFLAG) $(grbFPPFLAG) $(memFPPFLAG)
  # Extra flags
EXFLAGS =
###############################################################################
