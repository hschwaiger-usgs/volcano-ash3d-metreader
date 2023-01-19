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

# Script that converts the gfs grib files to netcdf using netcdf-java.
# This script is called from autorun_gfs.sh and takes three command-line arguments
#   convert_gfs.sh GFSres YYYYMMDD HR

# Select the netcdf version to write
#NCv=3
NCv=4

# Please edit these variables to match your system and location of netcdf-java
#JAVA="/usr/bin/java"
#NCJv="${HOME}/ncj/netcdfAll-4.5.jar"
rc=0
if [ $# -eq 0 ]
  then
  echo "No arguments supplied"
  echo "Usage: convert_gfs.sh Resolution YYYYMMDD FCpackage"
  echo "       where Resolution = 1p00, 0p50, or 0p25"
  echo "             YYYYMMDD   = year, month, day"
  echo "             FCpackage  = 00, 06, 12, or 18"
  exit
fi

echo "Looking for latest netcdfAll in ~/ncj/"
ls -1r ~/ncj/netcdfAll*.jar
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: Could not find netcdfAll in ~/ncj/ rc=$rc"
  echo "Please make sure to put netcdfAll-[].jar in ~/ncj or specify the path."
  echo "The latest version can be downloaded from"
  echo "  https://www.unidata.ucar.edu/downloads/netcdf-java/"
  exit 1
fi
NCJv=`ls -1r ~/ncj/netcdfAll*.jar | head -n 1`
echo "Found $NCJv"

echo "Looking for java"
which java
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: Could not find java in your path rc=$rc"
  exit 1
fi
JAVA=`which java`
echo "Found ${JAVA}"

WINDROOT="/data/WindFiles"
INSTALLDIR="/opt/USGS"

GFS=$1
yearmonthday=$2
FChour=$3

case ${GFS} in
 0p25)
  echo "GFS 0.25 degree"
  ;;
 0p50)
  echo "GFS 0.50 degree"
  ;;
 1p00)
  echo "GFS 1.00 degree"
  ;;
 *)
  echo "GFS product not recognized"
  echo "Valid values: 0p25, 0p50, 1p00"
  exit
esac

echo "------------------------------------------------------------"
echo "running convert_gfs.sh ${GFS} ${yearmonthday} ${FChour}"
echo `date`
echo "------------------------------------------------------------"

case ${GFS} in
 0p25)
  # GFS 0.25 degree
  HourMax=99
  HourStep=3
  #        gfs.t00z.pgrb2.0p25.f$000
  FilePre="gfs.t${FChour}z.pgrb2.0p25.f"
  iwf=22
  ;;
 0p50)
  # GFS 0.50 degree
  HourMax=198
  #HourMax=99
  HourStep=3
  #        gfs.t00z.pgrb2.0p50.f$000
  FilePre="gfs.t${FChour}z.pgrb2.0p50.f"
  iwf=20
  ;;
 1p00)
  # GFS 1.00 degree
  HourMax=384
  HourStep=3
  #        gfs.t00z.pgrb2.1p00.f$000
  FilePre="gfs.t${FChour}z.pgrb2.1p00.f"
  iwf=21
  ;;
 *)
  echo "GFS product not recognized"
  echo "Valid values: 0p25, 0p50, 1p00"
  exit
esac

rc=0
validlist="valid_files.txt"
GFSDATAHOME="${WINDROOT}/gfs"
#if [[ -d ${GFSDATAHOME} ]] ; then
#   echo "Error:  Download directory ${GFSDATAHOME} does not exist"
#   rc=$((rc + 1))
#   exit $rc
#fi

#name of directory containing current files
FC_day=gfs.${yearmonthday}${FChour}

#******************************************************************************
#START EXECUTING

#go to correct directory
echo "going to ${GFSDATAHOME}/${FC_day}"
cd ${GFSDATAHOME}/${FC_day}

#Convert to NetCDF
t=0
rm -f ${GFSDATAHOME}/${FC_day}/${validlist}
touch ${GFSDATAHOME}/${FC_day}/${validlist}
vcount=0
while [ "$t" -le ${HourMax} ]
do
  if [ "$t" -le 9 ]; then
      hour="00$t"
   elif [ "$t" -le 99 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  gfsfile=${FilePre}${hour}
  ncmlfile="gfs.t${FChour}z.f${hour}.ncml"
  netcdffile="${yearmonthday}${FChour}.f${hour}.nc"
  netcdf4file="${yearmonthday}${FChour}.f${hour}.nc4"
  if test -r ${gfsfile}
  then
     echo "making ${ncmlfile}"
     ${INSTALLDIR}/bin/makegfsncml ${gfsfile} ${ncmlfile}
     if [[ $? -ne 0 ]]; then
          exit 1
     fi
     echo "Converting ${gfsfile} to ${netcdffile}"
     # Note: nc4 is significantly smaller, but the direct conversion to nc4 from netcdf-java via the flag -netcdf4
     #       results in incompatible files.  Here, we use netcdf-java to product a temporary file, then convert
     #       to nc4 with nccopy
     echo "java -Xmx2048m -classpath ${NCJv} ucar.nc2.dataset.NetcdfDataset -in ${gfsfile} -out ${netcdffile} -IsLargeFile"
     if [ $NCv -eq 4 ]
     then
       ${JAVA} -Xmx2048m -classpath ${NCJv} ucar.nc2.dataset.NetcdfDataset -in ${gfsfile} -out tmp.nc -IsLargeFile
       nccopy -k 4 -d 5 tmp.nc ${netcdffile}
       rm tmp.nc
     else
       ${JAVA} -Xmx2048m -classpath ${NCJv} ucar.nc2.dataset.NetcdfDataset -in ${gfsfile} -out ${netcdffile} -IsLargeFile
     fi
     if [[ $? -ne 0 ]]; then
          exit 1
     fi
     # Check converted file for valid values
     echo "checking ${netcdffile} for corrupt values"
     ${INSTALLDIR}/bin/MetCheck ${iwf} 2 ${netcdffile}
     if [[ $? -eq 0 ]]; then
       cat MetCheck_log.txt >> ${GFSDATAHOME}/${FC_day}/${validlist}
       vcount=$((vcount+1))
     fi
   else
     echo "Warning: ${gfsfile} does not exist. Skipping."
   fi
   t=$((t+${HourStep}))
done

#Make sure the netcdf files all exist
echo "making sure all netcdf files exist"
t=0        # time index
numfiles=0 # file index
while [ "$t" -le ${HourMax} ]
do
  if [ "$t" -le 9 ]; then
      hour="00$t"
   elif [ "$t" -le 99 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  netcdffile="${yearmonthday}${FChour}.f${hour}.nc"
  if test -r ${netcdffile}
  then
     echo "${netcdffile} exists"
     t=$((t+${HourStep}))
     numfiles=$((numfiles+1))
   else
     echo "error: ${netcdffile} does not exist."
     exit 1
   fi
done

#set soft links in "latest" directory
echo "creating soft links in latest directory"
mkdir -p ${GFSDATAHOME}/latest
echo "rm ${GFSDATAHOME}/latest/*"
rm ${GFSDATAHOME}/latest/*
echo "rm ${GFSDATAHOME}/gfslist.txt"
rm ${GFSDATAHOME}/gfslist.txt
echo "4" > ${GFSDATAHOME}/gfslist.txt
echo "$numfiles" >> ${GFSDATAHOME}/gfslist.txt
t=0
while [ "$t" -le ${HourMax} ]
do
  if [ "$t" -le 9 ]; then
      hour="00$t"
   elif [ "$t" -le 99 ]; then
      hour="0$t"
   else
      hour="$t"
  fi
  netcdffile="${yearmonthday}${FChour}.f${hour}.nc"
  linkfile="latest.f${hour}.nc"
  echo "creating soft link for ${netcdffile}"
  ln -s ${GFSDATAHOME}/${FC_day}/${netcdffile} ${GFSDATAHOME}/latest/${linkfile}
  echo "latest/${linkfile}" >> ${GFSDATAHOME}/gfslist.txt
  t=$((t+3))
done

#cd ${GFSDATAHOME}
#if test -f "${INSTALLDIR}/bin/ncGFS4_2_pf"; then
#  PUFFDATAHOME="${WINDROOT}/puff/gfs/"
#  mkdir -p ${PUFFDATAHOME}
#  ${INSTALLDIR}/bin/ncGFS4_2_pf gfslist.txt
#  echo "mv Puff__GFS_______pf.nc ${PUFFDATAHOME}/${yearmonthday}${FChour}_gfs.nc"
#  mv Puff__GFS_______pf.nc ${yearmonthday}${FChour}_gfs.nc
#  mv ${yearmonthday}${FChour}_gfs.nc ${PUFFDATAHOME}/${yearmonthday}${FChour}_gfs.nc
#fi

echo "removing *.ncml, *.ncx2, and *.gbx9 files"
echo "rm ${GFSDATAHOME}/${FC_day}/*.ncml ${GFSDATAHOME}/${FC_day}/*.gbx9"
rm -f ${GFSDATAHOME}/${FC_day}/*.ncml ${GFSDATAHOME}/${FC_day}/*.gbx9 ${GFSDATAHOME}/${FC_day}/*.ncx2

echo "writing last_downloaded.txt"
echo ${yearmonthday}${FChour} > ${GFSDATAHOME}/last_downloaded.txt

echo "all done with windfiles"

echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
echo "finished convert_gfs.sh ${GFS} ${yearmonthday} ${FChour}"
echo `date`
echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
