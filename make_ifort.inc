###############################################################################
##########  Intel Fortran Compiler  #############################################

FCHOME = /opt/intel/oneapi/compiler/latest/linux/
FC = $(FCHOME)/bin/intel64/ifort
COMPINC = -I./ -I$(FCHOME)/include
COMPLIBS = -L./ -L$(FCHOME)/lib -L$(FCHOME)/compiler/lib/intel64_lin
LIBS = $(COMPLIBS) $(COMPINC)

# Debugging flags
ifeq ($(RUN), DEBUG)
    FFLAGS = -debug full -g3 -O0 -CB -init:snan,arrays -warn all -gen-interfaces -traceback -check all -check bounds -ftrapuv
endif
# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g2 -pg
endif
# Production run flags
ifeq ($(RUN), OPT)
    FFLAGS = -O3 -unroll -unroll-aggressive -ftz -inline-forceinline -finline-functions -ipo -ip
endif

  # Preprocessing flags
FPPFLAGS =  -fpp $(ncFPPFLAG) $(grbFPPFLAG) $(memFPPFLAG)
  # Extra flags
EXFLAGS =
###############################################################################
