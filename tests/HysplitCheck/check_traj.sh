#!/bin/bash

# This script pulls a Hysplit forecast for a volcano (1 = Akutan)
# The text data is parsed and a command is prepared for MetTraj_F
# Hysplit and MetTraj use different physics (Hysplit is much more
# sophisticated) so we really cannot easily numerically compare
# the results, but we can plot them with plot_trajHysplit.m
# and see how far off MetTraj is from Hysplit.

# Get today's Hysplit run for volc 1
infile="volc1_GFS_006.txt"
wget https://www.ready.noaa.gov/data/web/trajectory/volcanoes/${infile}

FC1=`sed -n 2p ${infile} | tr -s ' ' | cut -d' ' -f6`
FC2=`sed -n 3p ${infile} | tr -s ' ' | cut -d' ' -f6`

head -11 ${infile} | tail -7 > tmp.txt
YY=`head -1 tmp.txt | tr -s ' ' | cut -d' ' -f2`
MM=`head -1 tmp.txt | tr -s ' ' | cut -d' ' -f3`
DD=`head -1 tmp.txt | tr -s ' ' | cut -d' ' -f4`
HH=`head -1 tmp.txt | tr -s ' ' | cut -d' ' -f5`
LAT=`head -1 tmp.txt | tr -s ' ' | cut -d' ' -f6`
LON=`head -1 tmp.txt | tr -s ' ' | cut -d' ' -f7`
# Get the starting levels in m
LEV1=`sed -n 1p tmp.txt | tr -s ' ' | cut -d' ' -f8`
LEV2=`sed -n 2p tmp.txt | tr -s ' ' | cut -d' ' -f8`
LEV3=`sed -n 3p tmp.txt | tr -s ' ' | cut -d' ' -f8`
LEV4=`sed -n 4p tmp.txt | tr -s ' ' | cut -d' ' -f8`
LEV5=`sed -n 5p tmp.txt | tr -s ' ' | cut -d' ' -f8`
LEV6=`sed -n 6p tmp.txt | tr -s ' ' | cut -d' ' -f8`
LEV7=`sed -n 7p tmp.txt | tr -s ' ' | cut -d' ' -f8`

YYYY="20${YY}"   # converting year to 4-digits, assuming this century
HHH="${HH}.0"    # converting hours to floating-point
# Converting starting altitudes to km
LEV1km=$(echo "scale=3; $LEV1 / 1000.0" | bc -l)
LEV2km=$(echo "scale=3; $LEV2 / 1000.0" | bc -l)
LEV3km=$(echo "scale=3; $LEV3 / 1000.0" | bc -l)
LEV4km=$(echo "scale=3; $LEV4 / 1000.0" | bc -l)
LEV5km=$(echo "scale=3; $LEV5 / 1000.0" | bc -l)
LEV6km=$(echo "scale=3; $LEV6 / 1000.0" | bc -l)
LEV7km=$(echo "scale=3; $LEV7 / 1000.0" | bc -l)

# Reformat Hysplit data
dln=12  # start of data section of volc1_GFS_006.txt
nlv=7   # number of levels
nt=7    # number of time steps (including t=0)
# loop over time steps
rm -f htraj*.dat
touch htraj1.dat
touch htraj2.dat
touch htraj3.dat
touch htraj4.dat
touch htraj5.dat
touch htraj6.dat
touch htraj7.dat
for (( t=1;t<=$nt;t++))
do
  for (( l=1;l<=$nlv;l++))
  do
    ln=`echo "$dln + ($t-1)*$nlv + $l" | bc -l`
    lat=`sed -n ${ln}p ${infile} | tr -s ' ' | cut -d' ' -f11`
    lon=`sed -n ${ln}p ${infile} | tr -s ' ' | cut -d' ' -f12`
    echo "$lon $lat" >> htraj${l}.dat
  done
done

echo "/opt/USGS/bin/MetTraj_F ${LON} ${LAT} ${YYYY} ${MM} ${DD} ${HHH} 6.0 7 ${LEV1km} ${LEV2km} ${LEV3km} ${LEV4km} ${LEV5km} ${LEV6km} ${LEV7km} "
/opt/USGS/bin/MetTraj_F ${LON} ${LAT} ${YYYY} ${MM} ${DD} ${HHH} 6.0 7 ${LEV1km} ${LEV2km} ${LEV3km} ${LEV4km} ${LEV5km} ${LEV6km} ${LEV7km}

rm map_range_traj.txt tmp.txt
echo "Hysplit is using Forecast packages ${FC1} and ${FC2}"
echo "Make sure that MetTraj_F is using the correct windfile for comparison."
echo "To plot the results, run the octave script plot_trajHysplit.m"


