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

# Shell script that downloads hrrr data files (CONUS or AK) for the date supplied
# on the command line.
# This script is called from autorun_hrrr.sh and takes three command-line arguments
#   get_hrrr.sh HRRR YYYYMMDD HR

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

HRRR=$1
yearmonthday=$2
FChour=$3
SERVER="https://nomads.ncep.noaa.gov/pub/data/nccf/com/hrrr/prod"
WGETOPT="--no-check-certificate --tries=50"


echo "------------------------------------------------------------"
echo "running get_hrrr.sh script for $yearmonthday ${FChour}"
echo `date`
echo "------------------------------------------------------------"
t0=`date`

case ${HRRR} in
 1)
  # CONUS
  HourMax=18
  HourStep=1
  #        hrrr.t00z.wrfprsf00.grib2
  DirPre="conus"
  FilePre="hrrr.t${FChour}z.wrfprsf"
  FilePost=".grib2"
  ;;
 2)
  # AK
  HourMax=40
  HourStep=1
#        hrrr.t00z.wrfprsf00.ak.grib2
  DirPre="ak"
  FilePre="hrrr.t${FChour}z.wrfprsf"
  FilePost=".ak.grib2"
  ;;
 *)
  echo "HRRR product not recognized"
  echo "Valid values: 1 for CONUS or 2 for AK"
  exit
esac

HRRRDATAHOME="${WINDROOT}/hrrr/${DirPre}"
install -d ${HRRRDATAHOME}
if [[ $? -ne 0 ]] ; then
   echo "Error:  Download directory ${HRRRDATAHOME} cannot be"
   echo "        created or has insufficient write permissions."
   rc=$((rc + 1))
   exit $rc
fi

#name of directory containing current files
FC_day=${yearmonthday}_${FChour}

#******************************************************************************
#START EXECUTING

#go to correct directory
cd $HRRRDATAHOME
mkdir -p $FC_day
cd $FC_day

t=0
while [ "$t" -le ${HourMax} ]; do
  if [ "$t" -le 9 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  INFILE=${FilePre}${hour}${FilePost}
  fileURL=${SERVER}/hrrr.${yearmonthday}/${DirPre}/$INFILE
  echo "wget ${WGETOPT} ${fileURL}"
  time wget ${WGETOPT} ${fileURL}
  ${USGSROOT}/bin/gen_GRIB_index $INFILE
  ${USGSROOT}/bin/autorun_scripts/grib2nc.sh $INFILE
  t=$(($t+${HourStep}))
done

mkdir -p $HRRRDATAHOME/latest
cd $HRRRDATAHOME/latest
rm hrrr.*
ln -s ../$FC_day/* .

t1=`date`
echo "download start: $t0"
echo "download   end: $t1"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished get_hrrr.sh ${yearmonthday} ${FChour}"
echo `date`
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
