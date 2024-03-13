#!/bin/bash

sh clean.sh

# Looking for programs: bc and awk
rc=0
which bc > /dev/null
rc=$((rc + $?))
which awk > /dev/null
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Warning: Could not find bc or awk in your path"
  echo "         bc/awk are needed to verify accuracy of routines."
  echo "         Results will be compared with cmp"
  echo "         Values can be off by 0.001% using different compilation flags."
  HavBCawk=1
 else
  HavBCawk=0
fi

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

THRESH=0.001
if [ -z ${WINDROOT} ];then
 # Standard Linux location
 WINDROOT="/data/WindFiles"
 # Mac
 #WINDROOT="/opt/data/WindFiles"
fi
WindLabel="NCEP"
ln -s ${WINDROOT}/NCEP NCEP

sh clean.sh
echo "-----------------------------------------------------------"
echo "RUNNING NCEP TEST CASE 1: MetProbe"
echo "-----------------------------------------------------------"
rc=0
s=1
outdir="output${WindLabel}${s}"
../bin/MetProbe NCEP 1 1 -169.9468 52.8217 F 4 1 2 3 5 5 25 2 1980 5 17 0.0 > /dev/null 2>&1 
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: MetProbe returned error code"
  exit 1
fi
###  Run check on NWP_prof.dat and outputNCEP1/NWP_prof.dat
# This file has 5 data column and 17 rows
newdatafile="NWP_prof.dat"
olddatafile="output${WindLabel}${s}/NWP_prof.dat"
if [[ "$HavBCawk" -eq 0 ]] ; then
  for ir in {1..17}; do
   for ic in {1..5};  do
    new=`awk -v l="$ir" -v c="$ic" 'NR==l {printf ("%6.3f",$c)}' ${newdatafile}`
    old=`awk -v l="$ir" -v c="$ic" 'NR==l {printf ("%6.3f",$c)}' ${olddatafile}`
    err=$(echo "sqrt((($new - $old)/$old)^2)" | bc -l)
    st=`echo "$err > $THRESH" | bc`
    if [ $st -eq 1 ]; then
      # increment error code
      rc=$((rc + 1))
    fi
   done
  done
else
  cmp -s ${newdatafile} ${olddatafile}
  rc=$((rc + $?))
fi
if [[ $rc -eq 0 ]] ; then
  stat="PASS"
  printf " ---> ${GREEN}${stat}${NC}\n"
else
  stat="FAIL"
  printf " ---> ${RED}${stat}${NC}\n"
fi

sh clean.sh
echo "-----------------------------------------------------------"
echo "RUNNING NCEP TEST CASE 2: MetTraj_F"
echo "-----------------------------------------------------------"
rc=0
s=2
outdir="output${WindLabel}${s}"
cp ../examples/Traj_${WindLabel}.ctr .
../bin/MetTraj_F Traj_${WindLabel}.ctr  > /dev/null 2>&1
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: MetTraj_F returned error code"
  exit 1
fi

###  Run check on ftraj1.dat and outputSonde2/ftraj1.dat
# This file has 2 data column and 25 rows
newdatafile="ftraj1.dat"
olddatafile="output${WindLabel}${s}/ftraj1.dat"
if [[ "$HavBCawk" -eq 0 ]] ; then
  for ir in {1..25}; do
   for ic in {1..2};  do
    new=`awk -v l="$ir" -v c="$ic" 'NR==l {printf ("%6.3f",$c)}' ${newdatafile}`
    old=`awk -v l="$ir" -v c="$ic" 'NR==l {printf ("%6.3f",$c)}' ${olddatafile}`
    err=$(echo "sqrt((($new - $old)/$old)^2)" | bc -l)
    st=`echo "$err > $THRESH" | bc`
    if [ $st -eq 1 ]; then
      # increment error code
      rc=$((rc + 1))
    fi
   done
  done
else
  cmp -s ${newdatafile} ${olddatafile}
  rc=$((rc + $?))
fi

if [[ $rc -eq 0 ]] ; then
  stat="PASS"
  printf " ---> ${GREEN}${stat}${NC}\n"
else
  stat="FAIL"
  printf " ---> ${RED}${stat}${NC}\n"
fi

sh clean.sh
echo "-----------------------------------------------------------"
echo "RUNNING NCEP TEST CASE 3: MetRegrid 2d"
echo "-----------------------------------------------------------"
rc=0
s=3
ln -s ../examples/UpperAirSoundings .
outdir="output${WindLabel}${s}"
cp ../examples/Regrid_${WindLabel}.ctr .
../bin/MetRegrid Regrid_${WindLabel}.ctr  > /dev/null 2>&1
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: MetRegrid returned error code"
  exit 1
fi
###  Run check on out02_1_2.dat and outputSonde3/out02_1_2.dat
# This file has 7 data column and 22500 rows, but we only compare the sum of the last column
newdatafile="out02_1_2.dat"
olddatafile="output${WindLabel}${s}/out02_1_2.dat"
if [[ "$HavBCawk" -eq 0 ]] ; then
  ic=7
  new=`awk '{sum+=$7;} END{printf ("%6.3f",sum);}' ${newdatafile}`
  old=`awk '{sum+=$7;} END{printf ("%6.3f",sum);}' ${olddatafile}`
  err=$(echo "sqrt((($new - $old)/$old)^2)" | bc -l)
  st=`echo "$err > $THRESH" | bc`
  if [ $st -eq 1 ]; then
    # increment error code
    rc=$((rc + 1))
  fi
else
  cmp -s ${newdatafile} ${olddatafile}
  rc=$((rc + $?))
fi

if [[ $rc -eq 0 ]] ; then
  stat="PASS"
  printf " ---> ${GREEN}${stat}${NC}\n"
else
  stat="FAIL"
  printf " ---> ${RED}${stat}${NC}\n"
fi

rm -f NCEP

