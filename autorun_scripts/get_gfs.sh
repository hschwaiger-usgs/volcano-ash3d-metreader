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

# Shell script that downloads gfs data files (1, 0.5, or 0,25 deg.) for the date supplied
# on the command line.
# This script is called from autorun_gfs.sh and takes two command-line arguments
#   get_gfs.sh YYYYMMDD HR

# This is the location where metreader will be installed
INSTALLDIR="/opt/USGS"

# This is the location where the downloaded windfiles will be placed.
# Please edit this to suit your system.
if [ -z ${WINDROOT} ];then
 # Standard Linux location
 WINDROOT="/data/WindFiles"
 # Mac
 #WINDROOT="/opt/data/WindFiles"
fi
GFSDATAHOME="${WINDROOT}/gfs"

if [ $# -eq 0 ]
  then
  echo "No arguments supplied"
  echo "Usage: get_gfs.sh Resolution YYYYMMDD FCpackage"
  echo "       where Resolution = 1p00, 0p50, or 0p25"
  echo "             YYYYMMDD   = date"
  echo "             FCpackage  = 0, 6, 12, or 18"
  exit
fi

GFS=$1
yearmonthday=$2
FChour=$3

#SERVER="https://nomads.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"
SERVER="ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod"


echo "------------------------------------------------------------"
echo "running get_gfs.sh ${GFS} ${yearmonthday} ${FChour}"
echo `date`
echo "------------------------------------------------------------"
t0=`date`

case ${GFS} in
 0p25)
  # GFS 0.25 degree
  HourMax=99
  HourStep=3
  #        gfs.t00z.pgrb2.0p25.f$000
  FilePre="gfs.t${FChour}z.pgrb2.0p25.f"
  ;;
 0p50)
  # GFS 0.50 degree
  HourMax=198
  HourStep=3
  #        gfs.t00z.pgrb2.0p50.f$000
  FilePre="gfs.t${FChour}z.pgrb2.0p50.f"
  ;;
 1p00)
  # GFS 1.00 degree
  HourMax=384
  HourStep=3
  #        gfs.t00z.pgrb2.1p00.f$000
  FilePre="gfs.t${FChour}z.pgrb2.1p00.f"
  ;;
 *)
  echo "GFS product not recognized"
  echo "Valid values: 0p25, 0p50, 1p00"
  exit
esac


rc=0
GFSDATAHOME="${WINDROOT}/gfs"
install -d ${GFSDATAHOME}
if [[ $? -ne 0 ]] ; then
   echo "Error:  Download directory ${GFSDATAHOME} cannot be"
   echo "        created or has insufficient write permissions."
   rc=$((rc + 1))
   exit $rc
fi

#name of directory containing current files
FC_day=gfs.${yearmonthday}${FChour}

#******************************************************************************
#START EXECUTING

#go to correct directory
cd $GFSDATAHOME
mkdir -p $FC_day
cd $FC_day

t=0
while [ "$t" -le ${HourMax} ]; do
  if [ "$t" -le 9 ]; then
      hour="00$t"
   elif [ "$t" -le 99 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  INFILE=${FilePre}${hour}
  fileURL=${SERVER}/gfs.${yearmonthday}/${FChour}/atmos/$INFILE
  echo "wget ${fileURL}"
  time wget ${fileURL}
  ${INSTALLDIR}/bin/gen_GRIB_index $INFILE

  t=$(($t+${HourStep}))
done

echo "finished downloading wind files"
t1=`date`
echo "download start: $t0"
echo "download   end: $t1"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished get_gfs0.5deg.sh ${yearmonthday} ${FChour}"
echo `date`
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
