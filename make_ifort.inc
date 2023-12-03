###############################################################################
##########  Intel Fortran Compiler  #############################################
FCHOME = /opt/intel/oneapi/compiler/latest/linux/
FC = $(FCHOME)/bin/intel64/ifort
COMPINC = -I./ -I$(FCHOME)/include
COMPLIBS = -L./ -L$(FCHOME)/lib -L$(FCHOME)/compiler/lib/intel64_lin
LIBS = $(COMPLIBS) $(COMPINC)

# Debugging flags
ifeq ($(RUN), DEBUG)
    FFLAGS = -g2 -pg -warn all -check all -real-size 64 -check uninit -traceback -ftrapuv -debug all 
endif
# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g2 -pg
endif
# Production run flags
ifeq ($(RUN), OPT)
    FFLAGS = -O3 -ftz -w -ipo
endif

  # Preprocessing flags
FPPFLAGS =  -fpp $(ncFPPFLAG) $(grbFPPFLAG) $(memFPPFLAG)
  # Extra flags
EXFLAGS =
###############################################################################
