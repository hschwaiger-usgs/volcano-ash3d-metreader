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

# Shell script that downloads core data files for the date supplied on the command line.
# This script is takes four command-line arguments
#   get_core.sh YYYY MM DD
# This script uses the python download script provided by NOAA at:
#  https://ftp.cpc.ncep.noaa.gov/CORe/get_core/get_core.py
#
# Check environment variables WINDROOT and USGSROOT
#  WINDROOT = location where the downloaded windfiles will be placed.
#  USGSROOT = location where the MetReader tools and scripts were placed.
# Please edit these to suit your system or ensure WINDROOT/USGSROOT are set as environment
# variables in ${HOME}/.bash_profile or ${HOME}/.bashrc
if [ -z ${WINDROOT} ];then
 # default location
 WINDROOT="/data/WindFiles"
fi
if [ -z ${USGSROOT} ];then
 # default location
 USGSROOT="/opt/USGS"
fi
SCRIPTDIR="${USGSROOT}/bin/autorun_scripts"

if [ -z $3 ]
then
  echo "Error: Insufficient command-line arguments"
  echo "Usage:  get_core.sh YYYY MM DD"
  exit 1
else
  YYYY=$1
  MM=$2
  DD=$3
  yyyymmdd=${YYYY}${MM}${DD}
  echo "Downloading core files for date : ${yyyymmdd}"
fi

rc=0
COREDATAHOME="${WINDROOT}/CORe"
install -d ${COREDATAHOME}
if [[ $? -ne 0 ]]; then
   echo "Error:  Download directory ${COREDATAHOME} cannot be"
   echo "        created or has insufficient write permissions."
   rc=$((rc + 1))
   exit $rc
fi
echo "COREDATAHOME=$COREDATAHOME"

#if the directory for this year doesn't exist (e.g. it's Jan. 1), create it
echo "making sure the directory for year ${YYYY} exists"
if [ ! -r "${COREDATAHOME}/${YYYY}" ]
then
   echo "It doesnt.  Creating directory for year ${YYYY}"
   mkdir ${COREDATAHOME}/${YYYY}
else
   echo "Good.  It does."
fi

${SCRIPTDIR}/get_core.py pgb ${yyyymmdd}00 1 1 '(HGT|UGRD|VGRD|VVEL|TMP):(1000|925|850|800|750|700|600|500|400|300|250|200|150|100|70|50|30|20|10|5|2|1) mb' ${COREDATAHOME}/${YYYY}
${SCRIPTDIR}/get_core.py pgb ${yyyymmdd}03 1 1 '(HGT|UGRD|VGRD|VVEL|TMP):(1000|925|850|800|750|700|600|500|400|300|250|200|150|100|70|50|30|20|10|5|2|1) mb' ${COREDATAHOME}/${YYYY}
${SCRIPTDIR}/get_core.py pgb ${yyyymmdd}06 1 1 '(HGT|UGRD|VGRD|VVEL|TMP):(1000|925|850|800|750|700|600|500|400|300|250|200|150|100|70|50|30|20|10|5|2|1) mb' ${COREDATAHOME}/${YYYY}
${SCRIPTDIR}/get_core.py pgb ${yyyymmdd}09 1 1 '(HGT|UGRD|VGRD|VVEL|TMP):(1000|925|850|800|750|700|600|500|400|300|250|200|150|100|70|50|30|20|10|5|2|1) mb' ${COREDATAHOME}/${YYYY}
${SCRIPTDIR}/get_core.py pgb ${yyyymmdd}12 1 1 '(HGT|UGRD|VGRD|VVEL|TMP):(1000|925|850|800|750|700|600|500|400|300|250|200|150|100|70|50|30|20|10|5|2|1) mb' ${COREDATAHOME}/${YYYY}
${SCRIPTDIR}/get_core.py pgb ${yyyymmdd}15 1 1 '(HGT|UGRD|VGRD|VVEL|TMP):(1000|925|850|800|750|700|600|500|400|300|250|200|150|100|70|50|30|20|10|5|2|1) mb' ${COREDATAHOME}/${YYYY}
${SCRIPTDIR}/get_core.py pgb ${yyyymmdd}18 1 1 '(HGT|UGRD|VGRD|VVEL|TMP):(1000|925|850|800|750|700|600|500|400|300|250|200|150|100|70|50|30|20|10|5|2|1) mb' ${COREDATAHOME}/${YYYY}
${SCRIPTDIR}/get_core.py pgb ${yyyymmdd}21 1 1 '(HGT|UGRD|VGRD|VVEL|TMP):(1000|925|850|800|750|700|600|500|400|300|250|200|150|100|70|50|30|20|10|5|2|1) mb' ${COREDATAHOME}/${YYYY}

# Use this to convert to netcdf if you want
#/opt/USGS/bin/autorun_scripts/grib2nc.sh pgb.${yyyymmdd}00.grb 

