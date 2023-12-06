###############################################################################
##########  AMD Optimizing C/C++/Fortran Compiler (aocc) ######################

FCHOME = /opt/AMD/aocc-compiler-4.0.0
FC = $(FCHOME)/bin/flang
COMPINC = -I$(FCHOME)/include
COMPLIBS = -L./ -L$(FCHOME)/lib
LIBS = $(COMPLIBS) $(COMPINC)

# Debugging flags
ifeq ($(RUN), DEBUG)
    FFLAGS = -O0 -g3 -Wall -Mstandard -fPIC -v
endif
# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g2 -pg -fPIC
endif
# Production run flags
ifeq ($(RUN), OPT)
    #FFLAGS = -O3 -zopt -ftz -w -ipo
    FFLAGS = -Ofast -fPIC -mavx -fplugin-arg-dragonegg-llvm-codegen-optimize=3 -fplugin-arg-dragonegg-llvm-ir-optimize=3
endif

  # Preprocessing flags
FPPFLAGS =  -Mpreprocess $(ncFPPFLAG) $(grbFPPFLAG) $(memFPPFLAG)
  # Extra flags
EXFLAGS =
###############################################################################
