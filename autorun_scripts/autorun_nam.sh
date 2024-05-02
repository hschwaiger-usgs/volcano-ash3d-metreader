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

# Shell script that manages the download of the nam data files (AK 2.95 km or HI 2.5 km) for the
# current date.
# This script expects a command line argument indicating which forecast package to download.
#   autorun_nam.sh 091 0   for the AK HiRes 00 forecast package

# Please edit the line below to be consistant with the install directory specified in
# the makefile
INSTALLDIR="/opt/USGS"
ABRIDGED="0"            # Default is to get the whole file

if [ $# -eq 0 ]
  then
  echo "No arguments supplied"
  echo "Usage: autorun_nam.sh 091 0"
  exit
fi

NAM=$1
FC=$2

case ${NAM} in
 181)
  echo "Caribbean 0.108 deg"
  ;;
 196)
  echo "HI 2.5 km"
  ;;
 091)
  echo "AK 2.95"
  ABRIDGED="1"  # These files are too big so use a special script to pull
                # only select grib records
  ;;
 *)
  echo "NAM product not recognized"
  echo "Valid values: 091, 196"
  exit
esac

case ${FC} in
 0)
  FChour="00"
  FChourR="0.0"
  ;;
 6)
  FChour="06"
  FChourR="6.0"
  ;;
 12)
  FChour="12"
  FChourR="12.0"
  ;;
 18)
  FChour="18"
  FChourR="18.0"
  ;;
 *)
  echo "NAM forecast package not recognized"
  echo "Valid values: 0, 6, 12, 18"
  exit
esac

yearmonthday=`date -u +%Y%m%d`

echo "------------------------------------------------------------"
echo "running autorun_nam script : ${NAM} ${yearmonthday} ${FChour}"
echo "------------------------------------------------------------"

SCRIPTDIR="${INSTALLDIR}/bin/autorun_scripts"

#script that gets the wind files
if [ "$ABRIDGED" -eq "1" ]; then
 echo "  Calling ${SCRIPTDIR}/get_nam${NAM}.sh ${yearmonthday} ${FChour}"
 ${SCRIPTDIR}/get_nam${NAM}.sh ${yearmonthday} ${FChour}
else
 echo "  Calling ${SCRIPTDIR}/get_nam.sh ${NAM} ${yearmonthday} ${FChour}"
 ${SCRIPTDIR}/get_nam.sh ${NAM} ${yearmonthday} ${FChour}
fi
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished autorun_nam script"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
