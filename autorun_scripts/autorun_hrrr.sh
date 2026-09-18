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

# Shell script that manages the download of the hrrr data files (CONUS 3.0 km or AK 3.0 km) for the
# current date.
# This script expects a command line argument indicating which forecast package to download.
#   autorun_hrrr.sh 2 0   for the AK HiRes 00 forecast package
#
# or from a cron job:
#01 10 * * * /opt/USGS/bin/autorun_scripts/autorun_hrrr.sh 2 0    > /home/ash3d/cron_logs/hrrr200_log   2>&1
#
# Check environment variable USGSROOT
#  USGSROOT = location where the MetReader tools and scripts were placed.
# Please edit these to suit your system or ensure USGSROOT is set as environment
# variables in ${HOME}/.bash_profile or ${HOME}/.bashrc
if [ -z ${USGSROOT} ];then
 # default location
 USGSROOT="/opt/USGS"
fi

ABRIDGED="0"            # Default is to get the whole file

if [ $# -eq 0 ]
  then
  echo "No arguments supplied"
  echo "Usage: autorun_hrrr.sh 2 0"
  exit
fi

DOM=$1
FC=$2

case ${DOM} in
 1)
  echo "CONUS"
  ;;
 2)
  echo "AK"
  ;;
 *)
  echo "hrrr product not recognized"
  echo "Valid values: 1 (CONUS), 2 (AK)"
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

SCRIPTDIR="${USGSROOT}/bin/autorun_scripts"

#script that gets the wind files
echo "  Calling ${SCRIPTDIR}/get_hrrr.sh ${DOM} ${yearmonthday} ${FChour}"
${SCRIPTDIR}/get_hrrr.sh ${DOM} ${yearmonthday} ${FChour}
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished autorun_nam script"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
