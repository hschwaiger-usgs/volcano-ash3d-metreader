##############################################################################
#  Makefile for libmetreader.a
#
#    User-specified flags are in this top block
#
###############################################################################

#      This file is a component of the volcanic ash transport and dispersion model Ash3d,
#      written at the U.S. Geological Survey by Hans F. Schwaiger (hschwaiger@usgs.gov),
#      Larry G. Mastin (lgmastin@usgs.gov), and Roger P. Denlinger (roger@usgs.gov).

#      The model and its source code are products of the U.S. Federal Government and therefore
#      bear no copyright.  They may be copied, redistributed and freely incorporated 
#      into derivative products.  However as a matter of scientific courtesy we ask that
#      you credit the authors and cite published documentation of this model (below) when
#      publishing or distributing derivative products.

#      Schwaiger, H.F., Denlinger, R.P., and Mastin, L.G., 2012, Ash3d, a finite-
#         volume, conservative numerical model for ash transport and tephra deposition,
#         Journal of Geophysical Research, 117, B04204, doi:10.1029/2011JB008968. 

#      We make no guarantees, expressed or implied, as to the usefulness of the software
#      and its documentation for any purpose.  We assume no responsibility to provide
#      technical support to users of this software.

#  SYSTEM specifies which compiler to use
#    Current available options are:
#      gfortran , ifort
#    This variable cannot be left blank
#      
SYSTEM = gfortran
#SYSTEM = ifort
#
#  RUN specifies which collection of compilation flags that should be run
#    Current available options are:
#      DEBUG : includes debugging info and issues warnings
#      PROF  : includes profiling flags with some optimization
#      OPT   : includes optimizations flags for fastest runtime
#    This variable cannot be left blank
RUN=DEBUG
#RUN=PROF
#RUN=OPT
#
INSTALLDIR=/opt/USGS
#INSTALLDIR=$(HOME)/intel
#
# DATA FORMATS
#  For each data format you want to include in the library, set the corresponding
#  variable below to 'T'.  Set to 'F' any you do not want compiled or any unavailable
USENETCDF = T
USEGRIB   = T

# MEMORY
# If you need pointer arrays instead of allocatable arrays, set this to 'T'
USEPOINTERS = F

###############################################################################
#####  END OF USER SPECIFIED FLAGS  ###########################################
###############################################################################

FPPFLAGS = 
ifeq ($(USENETCDF), T)
 ncFPPFLAG = -DUSENETCDF
 ncOBJS = MetReader_NetCDF.o
 nclib = -lnetcdf -lnetcdff
else
 ncFPPFLAG =
 ncOBJS =
 nclib =
endif
ifeq ($(USEGRIB), T)
 grbFPPFLAG = -DUSEGRIB
 grbOBJS = MetReader_GRIB.o MetReader_GRIB_index.o
 grblib = -leccodes -leccodes_f90
else
 grbFPPFLAG =
 grbOBJS =
 grblib =
endif

ifeq ($(USEPOINTERS), T)
 memFPPFLAG = -DUSEPOINTERS
else
 memFPPFLAG =
endif

# location of HoursSince and projection
USGSLIBDIR = -L$(INSTALLDIR)/lib
USGSINC = -I$(INSTALLDIR)/include
USGSLIB = $(USGSINC) $(USGSLIBDIR) -lhourssince -lprojection

###############################################################################
###############################################################################

###############################################################################
##########  GNU Fortran Compiler  #############################################
ifeq ($(SYSTEM), gfortran)
    FCHOME=/usr
    FC = /usr/bin/gfortran
    #COMPINC = -I$(FCHOME)/local/include -I$(FCHOME)/include -I$(FCHOME)/lib64/gfortran/modules -I$(INSTALLDIR)/include
    #COMPLIBS = -L$(FCHOME)/local/lib -L$(FCHOME)/lib64 -L${INSTALLDIR}/lib 
    COMPINC = -I./ -I$(FCHOME)/include -I$(FCHOME)/lib64/gfortran/modules -I$(INSTALLDIR)/include
    COMPLIBS = -L./ -L$(FCHOME)/lib64 -L${INSTALLDIR}/lib

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
endif
###############################################################################
##########  Intel Fortran Compiler  #############################################
ifeq ($(SYSTEM), ifort)
    FCHOME = /opt/intel/oneapi/compiler/latest/linux/
    FC = $(FCHOME)/bin/intel64/ifort
    COMPINC = -I./ -I$(FCHOME)/include
    COMPLIBS = -L./ -L$(FCHOME)/lib
    LIBS = $(COMPLIBS) $(COMPINC)

# Debugging flags
ifeq ($(RUN), DEBUG)
    FFLAGS = -g2 -pg -warn all -check all -real-size 32 -check uninit -traceback -ftrapuv -debug all
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
endif
###############################################################################


LIB = libMetReader.a

ifeq ($(USEGRIB), T)
  GRIBTOOL = gen_GRIB_index
else
  GRIBTOOL =
endif

EXEC = \
 bin/MetRegrid \
 bin/MetSonde  \
 bin/MetTraj_F \
 bin/MetTraj_B \
 bin/MetCheck  \
 bin/probe_Met \
 bin/makegfsncml \
 bin/$(GRIBTOOL)

AUTOSCRIPTS = \
 autorun_scripts/autorun_gfs.sh \
 autorun_scripts/get_gfs.sh     \
 autorun_scripts/convert_gfs.sh \
 autorun_scripts/autorun_nam.sh \
 autorun_scripts/get_nam.sh     \
 autorun_scripts/autorun_NCEP_50YearReanalysis.sh \
 autorun_scripts/get_NCEP_50YearReanalysis.sh     \
 autorun_scripts/grib2nc.sh \
 autorun_scripts/prune_windfiles.sh \
 autorun_scripts/get_gmao.sh \
 autorun_scripts/probe_volc.sh
###############################################################################
# TARGETS
###############################################################################

all: ${lib} tools

lib: $(LIB)

tools: $(EXEC)

libMetReader.a: MetReader.F90 MetReader.o $(ncOBJS) $(grbOBJS) MetReader_Grids.o MetReader_ASCII.o makefile
	ar rcs libMetReader.a MetReader.o $(ncOBJS) $(grbOBJS) MetReader_Grids.o MetReader_ASCII.o

MetReader.o: MetReader.F90 makefile
	./get_version.sh
	$(FC) $(FPPFLAGS) $(EXFLAGS) -c MetReader.F90
MetReader_Grids.o: MetReader_Grids.f90 MetReader.o makefile
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(USGSLIB) -c MetReader_Grids.f90
MetReader_ASCII.o: MetReader_ASCII.f90 MetReader.o makefile
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(USGSLIB) -c MetReader_ASCII.f90

ifeq ($(USENETCDF), T)
MetReader_NetCDF.o: MetReader_NetCDF.F90 MetReader.o makefile
	$(FC) $(FPPFLAGS) $(FFLAGS) $(EXFLAGS) $(LIBS) $(nclib) $(USGSLIB) -c MetReader_NetCDF.F90
endif
ifeq ($(USEGRIB), T)
MetReader_GRIB_index.o: MetReader_GRIB_index.f90 makefile
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(grblib) $(USGSLIB) -c MetReader_GRIB_index.f90
MetReader_GRIB.o: MetReader_GRIB.f90 MetReader_GRIB_index.o MetReader.o makefile
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(grblib) -c MetReader_GRIB.f90
bin/gen_GRIB_index: gen_GRIB_index.f90 MetReader_GRIB_index.o makefile libMetReader.a
	mkdir -p bin
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(grblib) -c gen_GRIB_index.f90
	$(FC) $(FFLAGS) $(EXFLAGS) MetReader_GRIB_index.o gen_GRIB_index.o $(LIBS) $(grblib) -o bin/gen_GRIB_index
endif

bin/MetRegrid: tools/MetRegrid.f90 makefile libMetReader.a
	mkdir -p bin
	$(FC)  $(FFLAGS) $(EXFLAGS) tools/MetRegrid.f90 -o bin/MetRegrid -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB)

bin/MetSonde: tools/MetSonde.f90 makefile libMetReader.a
	mkdir -p bin
	$(FC) $(FFLAGS) $(EXFLAGS) -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB) -c tools/MetSonde.f90
	$(FC) $(FFLAGS) $(EXFLAGS) MetSonde.o  -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB) -o bin/MetSonde
bin/MetTraj_F: tools/MetTraj.F90 makefile libMetReader.a
	mkdir -p bin
	$(FC) $(FPPFLAGS) -DFORWARD  $(FFLAGS) $(EXFLAGS) tools/MetTraj.F90 -o bin/MetTraj_F -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB)
bin/MetTraj_B: tools/MetTraj.F90 makefile libMetReader.a
	mkdir -p bin
	$(FC) $(FPPFLAGS) -DBACKWARD $(FFLAGS) $(EXFLAGS) tools/MetTraj.F90 -o bin/MetTraj_B -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB)
bin/MetCheck: tools/MetCheck.f90 makefile libMetReader.a
	mkdir -p bin
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(nclib) $(grblib) -c tools/MetCheck.f90
	$(FC) $(FFLAGS) $(EXFLAGS) MetCheck.o -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB) -o bin/MetCheck
bin/probe_Met: tools/probe_Met.f90 makefile libMetReader.a
	mkdir -p bin
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(nclib) $(grblib) -c tools/probe_Met.f90
	$(FC) $(FFLAGS) $(EXFLAGS) probe_Met.o -L./ -lMetReader $(LIBS) $(nclib) $(grblib) $(USGSLIB) -o bin/probe_Met
bin/makegfsncml: tools/makegfsncml.f90 makefile
	mkdir -p bin
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(nclib) -c tools/makegfsncml.f90
	$(FC) $(FFLAGS) $(EXFLAGS) makegfsncml.o  $(LIBS) $(nclib) -o bin/makegfsncml

clean:
	rm -f *.o *__genmod.f90 *__genmod.mod
	rm -f *.mod
	rm -f lib*.a
	rm -f $(EXEC)

install:
	install -d $(INSTALLDIR)/lib/
	install -d $(INSTALLDIR)/include/
	install -d $(INSTALLDIR)/bin/
	install -d $(INSTALLDIR)/bin/autorun_scripts
	install -d $(INSTALLDIR)/share
	install -m 644 $(LIB) $(INSTALLDIR)/lib/
	install -m 644 *.mod $(INSTALLDIR)/include/
	install -m 755 $(EXEC) $(INSTALLDIR)/bin/
	install -m 755 $(AUTOSCRIPTS) $(INSTALLDIR)/bin/autorun_scripts/
	install -m 644 share/volc_NOVAC.txt $(INSTALLDIR)/share/volc_NOVAC.txt

uninstall:
	rm -f $(INSTALLDIR)/lib/$(LIB)
	rm -f $(INSTALLDIR)/include/metreader.mod
	rm -f $(INSTALLDIR)/bin/MetRegrid
	rm -f $(INSTALLDIR)/bin/MetSonde
	rm -f $(INSTALLDIR)/bin/MetTraj_F
	rm -f $(INSTALLDIR)/bin/MetTraj_B
	rm -f $(INSTALLDIR)/bin/MetCheck
	rm -f $(INSTALLDIR)/bin/makegfsncml
	rm -f $(INSTALLDIR)/bin/gen_GRIB_index
	rm -f $(INSTALLDIR)/bin/probe_Met
	rm -f $(INSTALLDIR)/bin/autorun_scripts/autorun_gfs.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/autorun_nam.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/autorun_NCEP_50YearReanalysis.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/get_gfs.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/get_nam.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/get_NCEP_50YearReanalysis.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/convert_gfs.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/grib2nc.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/prune_windfiles.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/get_gmao.sh
	rm -f $(INSTALLDIR)/bin/autorun_scripts/probe_volc.sh
	rm -f $(INSTALLDIR)/share/volc_NOVAC.txt

