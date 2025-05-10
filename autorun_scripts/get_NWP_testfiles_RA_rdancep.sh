#!/bin/bash
  
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

# Shell script that downloads NASA MERRA-2 Reanalysis data files for the date supplied
# on the command line.

YYYY=$1
MM=$2
DD=$3

if [ $# -eq 0 ]
  then
  echo "No arguments supplied"
  echo "Usage: get_NWP_testfiles_RA_nasa.sh YYYY MM DD"
  exit
fi

# You can override the date here:
#YYYY=2018
#MM=06
#DD=20

YYYYMMDD=${YYYY}${MM}${DD}
# Where to store the downloads
DATAHOME="/data/NWP_testfiles/Reanalysis/NetCDF"
mkdir -p ${DATAHOME}/${YYYYMMDD}
cd ${DATAHOME}/${YYYYMMDD}

#---------------------------------------
#  ERA5
UCARSERVER="https://data-osdf.rda.ucar.edu/ncar/rda"
ERA5_remote="${UCARSERVER}/d633000/e5.oper.an.pl/${YYYY}/${MM}"
# This product is available after about 10 weeks from present
# ~ 211 Mb for 1-day of data
# iwf = 
wget ${ERA5_remote}/e5.oper.an.pl.128_248_cc.ll025sc.${YYYYMMDD}00_${YYYYMMDD}23.nc

#---------------------------------------
#  JRA-3Q
UCARSERVER="https://data-osdf.rda.ucar.edu/ncar/rda"
JRA_remote="${UCARSERVER}/d640000/anl_p/${YYYY}/${MM}"
# This product is available after about 2 weeks from present
# ~ 300 Mb for 5-days of data (need to sort out day indexing)
# iwf = 
wget ${JRA_remote}/jra3q.anl_p.0_3_5.hgt-pres-an-gauss.${YYYY}${MM}2600_${YYYY}${MM}3018.nc

#---------------------------------------
#  NARR
UCARSERVER="https://data-osdf.rda.ucar.edu/ncar/rda"
NARR_remote="${UCARSERVER}/d608000/NARRraw/${YYYY}${MM}"
# This product is available after about 2 weeks from present
# ~ 445 Mb for 1-day of data
# iwf = 3
wget ${NARR_remote}/merged_AWIP32.${YYYYMMDD}.tar



