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

# Parsing command-line arguments
#  first/second pass , rundirectory
echo "------------------------------------------------------------"
echo "running GMT_plot_traj.sh with parameter:"
if [ "$#" -eq 1 ]; then
 if [ $1 -eq 0 ]; then
   echo " Foreward trajectories requested"
   FORE=0
   pref='f'
 elif [ $1 -eq 1 ]; then
   echo " Backward trajectories requested"
   FORE=1
   pref='b'
 else
   echo " Argument not recognized."
   echo " Expected 0 for foreward trajectories"
   echo "          1 for backward trajectories"
   echo " Assuming foreward"
   FORE=0
   pref='f'
 fi
else
  echo "No command-line argument detected; assuming foreward trajectories."
  FORE=0
  pref='f'
fi

if [ "$#" -eq 2 ]; then
  echo "Second command line argument detected: setting run directory"
  RUNHOME=$2
 else
  echo "No second command line argument detected, using cwd"
  RUNHOME=`pwd`
fi
cd ${RUNHOME}
echo `date`
echo "------------------------------------------------------------"
rc=0                                             # error message accumulator
CLEANFILES="T"

# We need to know if we must prefix all gmt commands with 'gmt', as required by version 5/6
GMTv=5
type gmt >/dev/null 2>&1 || { echo >&2 "Command 'gmt' not found.  Assuming GMTv4."; GMTv=4;}
if [ $GMTv -eq 4 ] ; then
    echo "GMT 4 is no longer supported."
    echo "Please update to GMT 5 or 6"
    exit 1
 else
    GMTv=`gmt --version | cut -c1`
fi
GMTpen=("-" "-" "-" "-" "/" ",")
echo "GMT version = ${GMTv}"

USGSROOT="/opt/USGS"
ASH3DROOT="${USGSROOT}/Ash3d"

ASH3DBINDIR="${ASH3DROOT}/bin"
ASH3DSCRIPTDIR="${ASH3DROOT}/bin/scripts"
ASH3DSHARE="$ASH3DROOT/share"
ASH3DSHARE_PP="${ASH3DSHARE}/post_proc"

if test -r world_cities.txt ; then
    echo "Found file world_cities.txt"
  else
    ln -s ${ASH3DSHARE_PP}/world_cities.txt .
fi


# This is the initial run before the full Ash3d run
SUB=0
LLLON=`cat map_range_traj.txt | awk '{print $1}'`
LLLAT=`cat map_range_traj.txt | awk '{print $3}'`
URLON=`cat map_range_traj.txt | awk '{print $2}'`
URLAT=`cat map_range_traj.txt | awk '{print $4}'`
DLON=`echo "$URLON-$LLLON" | bc -l`
DLAT=`echo "$URLAT-$LLLAT" | bc -l`
# Now we need to adjust the limits so that the map has the approximately correct aspect ratio
dum=`echo "$DLAT * 2.0" | bc -l`
test1=`echo "$DLON < $dum" | bc -l`
echo "test1 = $test1"
#if (( $DLON < $dum )); then
if [ $test1 -eq 1 ]; then
  echo "Resetting DLON"
  DLON=$dum
  URLON=`echo "$LLLON+$DLON" | bc -l`
fi

echo "LLLON=$LLLON, LLLAT=$LLLAT, DLON=$DLON, DLAT=$DLAT"
echo "URLON=$URLON, URLAT=$URLAT"

VCLON=`head -n1 ${pref}traj1.dat | awk '{print $1}'`
VCLAT=`head -n1 ${pref}traj1.dat | awk '{print $2}'`

###############################################################################
##  Now make the maps

t=0
gmt gmtset PROJ_ELLIPSOID Sphere

AREA="-R$LLLON/$URLON/$LLLAT/$URLAT"
DLON_INT="$(echo $DLON | sed 's/\.[0-9]*//')"  #convert DLON to an integer
if [ $DLON_INT -le 5 ] ; then
   BASE="-Ba1/a1"                  # label every 5 degress lat/lon
   DETAIL="-Dh"                        # high resolution coastlines (-Dc=crude)
 elif [ $DLON_INT -le 10 ] ; then
   BASE="-Ba2/a2"                  # label every 5 degress lat/lon
   DETAIL="-Dh"                        # high resolution coastlines (-Dc=crude)
 elif [ $DLON_INT -le 20 ] ; then
   BASE="-Ba5/a5"                  # label every 5 degress lat/lon
   DETAIL="-Dh"                        # high resolution coastlines (-Dc=crude)
 else
   BASE="-Ba10/a10"                    #label every 10 degrees lat/lon
   DETAIL="-Dl"                        # low resolution coastlines (-Dc=crude)
fi
PROJ="-JM${VCLON}/${VCLAT}/20"
COAST="-G220/220/220 -W"            # RGB values for land areas (220/220/220=light gray)
BOUNDARIES="-Na"                    # -N=draw political boundaries, a=all national, Am. state & marine b.

#############################################################################
### Plot the base map
echo "gmt pscoast $AREA $PROJ $BASE $DETAIL $COAST $BOUNDARIES -K  > temp.ps"
gmt pscoast $AREA $PROJ $BASE $DETAIL $COAST $BOUNDARIES -K  > temp.ps

# This is the Hysplit-like trajectory plot using the same basemap
# Trajectories currently (2017-04-11) are:
#  5000 ft (1.5240 km) Red       (255/0/0)
# 10000 ft (3.0480 km) Blue      (0/0/255)
# 15000 ft (4.5720 km) Green     (0/255/0)    NOTE: This is not plotted currently
# 20000 ft (6.0960 km) Cyan      (0/255/255)
# 30000 ft (9.1440 km) Magenta   (255/0/255)
# 40000 ft (12.192 km) Yellow    (255/255/0)
# 50000 ft (15.240 km) Blue-grey (51/153/204)
gmt psxy ${pref}traj1.dat   $AREA $PROJ -P -K -O -W4${GMTpen[GMTv]}255/0/0    -V >> temp.ps
gmt psxy ${pref}traj2.dat   $AREA $PROJ -P -K -O -W4${GMTpen[GMTv]}0/0/255    -V >> temp.ps
gmt psxy ${pref}traj3.dat   $AREA $PROJ -P -K -O -W4${GMTpen[GMTv]}0/255/0    -V >> temp.ps
gmt psxy ${pref}traj4.dat   $AREA $PROJ -P -K -O -W4${GMTpen[GMTv]}0/255/255  -V >> temp.ps
gmt psxy ${pref}traj5.dat   $AREA $PROJ -P -K -O -W4${GMTpen[GMTv]}255/0/255  -V >> temp.ps
gmt psxy ${pref}traj6.dat   $AREA $PROJ -P -K -O -W4${GMTpen[GMTv]}255/255/0  -V >> temp.ps
gmt psxy ${pref}traj7.dat   $AREA $PROJ -P -K -O -W4${GMTpen[GMTv]}51/153/204 -V >> temp.ps

awk '{print $1, $2, 1.0}' ${pref}traj1.dat | gmt psxy $AREA $PROJ $BASE -Sc0.10i -W1${GMTpen[GMTv]}0/0/0 -G255/0/0    -O -K >> temp.ps
awk '{print $1, $2, 1.0}' ${pref}traj2.dat | gmt psxy $AREA $PROJ $BASE -Sc0.10i -W1${GMTpen[GMTv]}0/0/0 -G0/0/255    -O -K >> temp.ps
awk '{print $1, $2, 1.0}' ${pref}traj3.dat | gmt psxy $AREA $PROJ $BASE -Sc0.10i -W1${GMTpen[GMTv]}0/0/0 -G0/255/0    -O -K >> temp.ps
awk '{print $1, $2, 1.0}' ${pref}traj4.dat | gmt psxy $AREA $PROJ $BASE -Sc0.10i -W1${GMTpen[GMTv]}0/0/0 -G0/255/255  -O -K >> temp.ps
awk '{print $1, $2, 1.0}' ${pref}traj5.dat | gmt psxy $AREA $PROJ $BASE -Sc0.10i -W1${GMTpen[GMTv]}0/0/0 -G255/0/255  -O -K >> temp.ps
awk '{print $1, $2, 1.0}' ${pref}traj6.dat | gmt psxy $AREA $PROJ $BASE -Sc0.10i -W1${GMTpen[GMTv]}0/0/0 -G255/255/0  -O -K >> temp.ps
awk '{print $1, $2, 1.0}' ${pref}traj7.dat | gmt psxy $AREA $PROJ $BASE -Sc0.10i -W1${GMTpen[GMTv]}0/0/0 -G51/153/204 -O -K >> temp.ps
echo "Finished plotting trajectory data"

# Last gmt command is to plot the volcano and close out the ps file
echo $VCLON $VCLAT '1.0' | gmt psxy $AREA $PROJ -St0.1i -Gblack -Wthinnest -O >> temp.ps

#  Convert to gif
if [ $GMTv -eq 5 ] ; then
    gmt psconvert temp.ps -A -Tg
    convert -rotate 90 temp.png -resize 630x500 -alpha off temp.gif
  else
    gmt psconvert temp.ps -A -Tg
    convert temp.png -resize 630x500 -alpha off temp.gif
fi

mv temp.gif trajectory_${SUB}.gif

# Clean up temporary files
#if [ "$CLEANFILES" == "T" ]; then
#   rm -f *.grd *.lev
#   rm -f caption*.txt cities.xy map_range*txt legend_positions*txt
#   rm -f temp.*
#   rm -f gmt.conf gmt.history
#   rm -f world_cities.txt
#   rm -f legend*png
#fi

echo "exiting GMT_plot_traj.sh with status $rc"
exit $rc

