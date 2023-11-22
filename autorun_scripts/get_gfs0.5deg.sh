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

# Shell script that downloads gfs 0.5degree data files for the current date.
# This script is called from autorun_gfs0.5deg.sh and takes two command-line arguments
#   get_gfs0.5deg.sh YYYYMMDD HR

# This is the location where the downloaded windfiles will be placed.
# Please edit this to suit your system.
if [ -z ${WINDROOT} ];then
 # Standard Linux location
 WINDROOT="/data/WindFiles"
 # Mac
 #WINDROOT="/opt/data/WindFiles"
fi

yearmonthday=$1
FChour=$2

echo "------------------------------------------------------------"
echo "running get_gfs0.5deg.sh ${yearmonthday} ${FChour}"
echo `date`
echo "------------------------------------------------------------"
t0=`date`

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
mkdir $FC_day
cd $FC_day

t=0
while [ "$t" -le 99 ]; do
  if [ "$t" -le 9 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  filename="ftp://ftp.ncep.noaa.gov/pub/data/nccf/com/gfs/prod/gfs.${yearmonthday}/${FChour}/gfs.t${FChour}z.pgrb2.0p50.f${hour}"
  echo "wget ${filename}"
  time wget ${filename}
  t=$(($t+3))
done

echo "finished downloading wind files"
t1=`date`
echo "download start: $t0"
echo "download   end: $t1"
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished get_gfs0.5deg.sh ${yearmonthday} ${FChour}"
echo `date`
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
