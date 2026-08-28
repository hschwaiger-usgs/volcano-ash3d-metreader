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

# Shell script that manages the download of the CORe data files for the
# current date. This script uses the python download script provided by NOAA at:
#  https://ftp.cpc.ncep.noaa.gov/CORe/get_core/get_core.py
#
# or from a cron job:
#01 23 * * * /opt/USGS/bin/autorun_scripts/autorun_core.sh  > /home/ash3d/cron_logs/core_log  2>&1
#
# Check environment variables WINDROOT and USGSROOT
#  WINDROOT = location where the downloaded windfiles will placed.
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

echo "------------------------------------------------------------"
echo "running autorun_core.sh"
echo `date`
echo "------------------------------------------------------------"

SCRIPTDIR="${USGSROOT}/bin/autorun_scripts"

starttime=`date`                  #record when we're starting the download
echo "starting autorun_core.sh at ${starttime}"

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

# Get today's date
y=`date +%Y`
monthnow=`date +%m`
daynow=`date +%d`
# Now get date for the-day-before-yesterday since CORe seems to have a 2-day latency
YYYY=`date --date="${y}/$monthnow/${daynow} - 2 days" +%Y`
MM=`date --date="${y}/$monthnow/${daynow} - 2 days" +%m`
DD=`date --date="${y}/$monthnow/${daynow} - 2 days" +%d`
echo "year=$YYYY"

#if the directory for this year doesn't exist (e.g. it's Jan. 1), create it
echo "making sure the directory for year ${YYYY} exists"
if [ ! -r "${COREDATAHOME}/${YYYY}" ]
then
   echo "It doesnt.  Creating directory for year ${YYYY}"
   mkdir ${COREDATAHOME}/${YYYY}
else
   echo "Good.  It does."
fi

yyyymmdd=${YYYY}${MM}${DD}

${SCRIPTDIR}/get_core.sh ${YYYY} ${MM} ${DD}

# Log the last date downloaded
echo ${yyyymmdd} > ${COREDATAHOME}/last_downloaded.txt   #write date of last download to text file

endtime=`date`

echo "download started at $starttime"
echo "download ended at $endtime"

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished autorun_core.sh"
echo `date`
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"

