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

NASASERVER="https://goldsmr5.gesdisc.eosdis.nasa.gov/data/MERRA2/M2T3NVMST.5.12.4"
NASA_remote="${NASASERVER}/${YYYY}/${MM}"
#WGETOPT="--no-check-certificate --tries=50"
WGETOPT="--load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies"
# This product is available after about 5 weeks from present
# ~ 2.1 Gb for 1-day of data
# iwf = 40
wget ${WGETOPT} ${NASA_remote}/MERRA2_400.tavg3_3d_mst_Nv.${YYYY}${MM}${DD}.nc4

