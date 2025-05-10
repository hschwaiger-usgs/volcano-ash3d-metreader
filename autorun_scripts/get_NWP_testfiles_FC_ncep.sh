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

# Shell script that downloads gfs/nam/ecmwf forecast data files for the date supplied
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
DATAHOME="/data/NWP_testfiles/Forecasts/GRIB"
mkdir -p ${DATAHOME}/${YYYYMMDD}
cd ${DATAHOME}/${YYYYMMDD}

FC=00    # which forecast package
FCend=6  # number of forecast hours to download (starting at 0)
         # This only replaces one character so this shouldn't be >9

YYYYMMDDHH=${YYYY}${MM}${DD}${FC}

# This site hold the past 10 days or so of forecast packages
NCEPSERVER="ftp://ftp.ncep.noaa.gov/pub/data/nccf/com"
#NCEPSERVER="https://nomads.ncep.noaa.gov/pub/data/nccf/com"
ECMWFSERVER="https://data.ecmwf.int/forecasts"

GFSROOT_remote="${NCEPSERVER}/gfs/prod/gfs.${YYYYMMDD}/${FC}/atmos/"
GDASROOT_remote="${NCEPSERVER}/gfs/prod/gdas.${YYYYMMDD}/${FC}/atmos/"
NAMROOT_remote="${NCEPSERVER}/nam/prod/nam.${YYYYMMDD}"
ECMWF_remote="${ECMWFSERVER}/${YYYYMMDD}/${FC}z/ifs/0p25/oper"

# Downloading all forecast files with 1-hour steps
FCincr=1
t=0
while [ "$t" -le "$FCend" ]; do
  wget ${GFSROOT_remote}/gfs.t${FC}z.pgrb2.0p25.f00$t                  # ~520 Mb : GFS  0.25   iwf = 22
  wget ${GFSROOT_remote}/gfs.t${FC}z.pgrb2b.0p25.f00$t                 # ~230 Mb : GFSb 0.25   iwf = 22
  wget ${GDASROOT_remote}/gdas.t${FC}z.pgrb2.0p25.f00${t}              # ~490 Mb : GDAS 0.25   iwf = 22

  wget ${NAMROOT_remote}/nam.t${FC}z.alaskanest.hiresf0${t}.tm00.grib2 # ~800 Mb : Grid 91 AK 2.976
  wget ${NAMROOT_remote}/nam.t${FC}z.awip3d0${t}.tm00.grib2            #  ~11 Mb : Grid 212 CONUS (40.6 km)
  wget ${NAMROOT_remote}/nam.t${FC}z.awip120${t}.tm00.grib2            #  ~30 Mb : Grid ?
  wget ${NAMROOT_remote}/nam.t${FC}z.awip320${t}.tm00.grib2            #  ~50 Mb : Grid 221 N.Amer 32.5
  wget ${NAMROOT_remote}/nam.t${FC}z.awphys0${t}.tm00.grib2            #  ~57 Mb : Grid 218 CONUS (12.2 km)
  wget ${NAMROOT_remote}/nam.t${FC}z.conusnest.hiresf0${t}.tm00.grib2  # ~930 Mb : Grid 227 CONUS (5.1 km)
  wget ${NAMROOT_remote}/nam.t${FC}z.firewxnest.hiresf0${t}.tm00.grib2 # ~145 Mb :
  wget ${NAMROOT_remote}/nam.t${FC}z.grbgrd0${t}.tm00.grib2            #   ~9 Mb : Grid 104 N.Amer 90.7
  wget ${NAMROOT_remote}/nam.t${FC}z.hawaiinest.hiresf0${t}.tm00.grib2 #  ~28 Mb : Grid 196 HI (2.5 km)
  wget ${NAMROOT_remote}/nam.t${FC}z.priconest.hiresf0${t}.tm00.grib2  #  ~70 Mb : Grid 194 Puerto Rico (2.5 km)

  t=$(($t+${FCincr}))
done

FCincr=3
t=0
while [ "$t" -le "$FCend" ]; do
  wget ${ECMWF_remote}/${YYYYMMDD}000000-${t}h-oper-fc.grib2    # ~130 Mb : ECMWF 0.25 iwf = 34
  wget ${GFSROOT_remote}/gfs.t${FC}z.pgrb2.0p50.f00${t}        # ~155 Mb : GFS 0.50   iwf = 20
  wget ${GFSROOT_remote}/gfs.t${FC}z.pgrb2.1p00.f00${t}        #  ~43 Mb : GFS 1.00   iwf = 21
  wget ${GFSROOT_remote}/gfs.t${FC}z.pgrb2full.0p50.f00${t}    # ~220 Mb : GFS 0.50   iwf = 20

  wget ${NAMROOT_remote}/nam.t${FC}z.afwaca0${t}.tm00.grib2    #  ~40 Mb : Grid 181 Caribbean (0.108 deg)
  wget ${NAMROOT_remote}/nam.t${FC}z.afwahi0${t}.tm00.grib2    #  ~23 mb : Grid 182 HI (0.108 deg)
  wget ${NAMROOT_remote}/nam.t${FC}z.awak3d0${t}.tm00.grib2    #  ~55 Mb : Grid 242 AK 11.25
  wget ${NAMROOT_remote}/nam.t${FC}z.awip200${t}.tm00.grib2    #   ~5 Mb : 
  wget ${NAMROOT_remote}/nam.t${FC}z.awipak0${t}.tm00.grib2    #   ~7 Mb : Grid 216 AK 45.0
  wget ${NAMROOT_remote}/nam.t${FC}z.awiphi0${t}.tm00.grib2    #   ~6 Mb : Grid 243 E.N.Pac (0.4 deg)
  wget ${NAMROOT_remote}/nam.t${FC}z.awp2420${t}.tm00.grib2    #  ~22 Mb : 
  wget ${NAMROOT_remote}/nam.t${FC}z.bgrdsf0${t}.tm00.grib2    #  ~47 Mb : 
  t=$(($t+${FCincr}))
done

FCincr=6
t=0
while [ "$t" -le "$FCend" ]; do
  wget ${NAMROOT_remote}/nam.t${FC}z.awp2110${t}.tm00.grib2    # ~1.2 Mb : Grid 211 CONUS (81.3 km)
  wget ${NAMROOT_remote}/nam.t${FC}z.bgrd3d0${t}.tm00.grib2    # ~482 Mb : 
  t=$(($t+${FCincr}))
done


##########################################################
#
#NAM Grid 005  No longer available
