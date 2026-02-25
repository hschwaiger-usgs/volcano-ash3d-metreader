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

# Shell script that downloads NASA forecast data files for the date supplied
# on the command line.

YYYY=$1
MM=$2
DD=$3

if [ $# -eq 0 ]
  then
  echo "No arguments supplied"
  echo "Usage: get_NWP_testfiles_FC.sh YYYY MM DD"
  exit
fi

# You can override the date here:
#YYYY=2018
#MM=06
#DD=20

YYYYMMDD=${YYYY}${MM}${DD}
# Where to store the downloads
DATAHOME="/data/NWP_testfiles/Forecasts/NetCDF"
mkdir -p ${DATAHOME}/${YYYYMMDD}
cd ${DATAHOME}/${YYYYMMDD}

FC=00    # which forecast package
FCend=6  # number of forecast hours to download (starting at 0)
         # This only replaces one character so this shouldn't be >9

NASASERVER="https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/forecast"
NASAROOT_remote="${NASASERVER}/Y${YYYY}/M${MM}/D${DD}/H00"
WGETOPT="--no-check-certificate --tries=50"
# Downloading all forecast files with 3-hour steps
FCincr=3
t=0
while [ "$t" -le "$FCend" ]; do

  wget ${WGETOPT} ${NASAROOT_remote}/GEOS.fp.fcst.inst3_3d_asm_Cp.${YYYYMMDD}_00+${YYYYMMDD}_0${t}00.V01.nc4   # ~181 Mb : NASA  0.625 x 0.5   iwf = 40
  wget ${WGETOPT} ${NASAROOT_remote}/GEOS.fp.fcst.inst3_3d_asm_Np.${YYYYMMDD}_00+${YYYYMMDD}_0${t}00.V01.nc4   # ~710 Mb : NASA  0.3125 x 0.25 iwf = 41

  t=$(($t+${FCincr}))
done

