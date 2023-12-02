###############################################################################
##########  GNU Fortran Compiler  #############################################
FCHOME=/usr
FC = $(FCHOME)/bin/gfortran
#COMPINC = -I$(FCHOME)/local/include -I$(FCHOME)/include -I$(FCHOME)/lib64/gfortran/modules -I$(INSTALLDIR)/include
#COMPLIBS = -L$(FCHOME)/local/lib -L$(FCHOME)/lib64 -L${INSTALLDIR}/lib 
COMPINC = -I./ -I$(FCHOME)/include -I$(FCHOME)/lib64/gfortran/modules -I$(INSTALLDIR)/include
#  -I$(FCHOME)/lib/x86_64-linux-gnu/fortran/gfortran-mod-15
COMPLIBS = -L./ -L$(FCHOME)/lib -L$(FCHOME)/lib64 -L${INSTALLDIR}/lib

LIBS = $(COMPLIBS) $(COMPINC)
# -lefence 

# Debugging flags
ifeq ($(RUN), DEBUG)
    FFLAGS = -O0 -g3 -Wall -Wextra -fimplicit-none  -Wall  -Wline-truncation  -Wcharacter-truncation  -Wsurprising  -Waliasing  -Wimplicit-interface  -Wunused-parameter  -fwhole-file  -fcheck=all  -std=f2008  -pedantic  -fbacktrace -Wunderflow -ffpe-trap=invalid,zero,overflow -fdefault-real-8
endif
# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g -pg -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -fdefault-real-8
endif
# Production run flags
ifeq ($(RUN), OPT)
    FFLAGS = -O3 -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -fdefault-real-8
endif
  # Preprocessing flags
FPPFLAGS = -x f95-cpp-input $(ncFPPFLAG) $(grbFPPFLAG) $(memFPPFLAG)
  # Extra flags
EXFLAGS =
###############################################################################

