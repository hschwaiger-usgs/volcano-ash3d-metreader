#!/bin/bash

sh clean.sh

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

WindLabel="Sonde"

sh clean.sh
echo "-----------------------------------------------------------"
echo "RUNNING SONDE TEST CASE 1: MetProbe"
echo "-----------------------------------------------------------"
s=0
outdir="output${WindLabel}${s}"
cp ../examples/UpperAirSoundings/UIL_1980051700_raw.dat .
../bin/MetProbe UIL_1980051700_raw.dat 1 1 -169.9468 52.8217 F 4 1 2 3 5 1 2 1 1980 5 17 0.0 > /dev/null 2>&1
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: MetProbe returned error code"
  exit 1
fi
###  Run check on NWP_prof.dat and outputSonde1/NWP_prof.dat
cmp -s NWP_prof.dat output${WindLabel}1/NWP_prof.dat
if [[ $? == 0 ]] ; then
  stat="PASS"
  printf " ---> ${GREEN}${stat}${NC}\n"
else
  stat="FAIL"
  printf " ---> ${RED}${stat}${NC}\n"
fi

sh clean.sh
echo "-----------------------------------------------------------"
echo "RUNNING SONDE TEST CASE 2: MetTraj_F"
echo "-----------------------------------------------------------"
s=1
ln -s ../examples/UpperAirSoundings .
outdir="output${WindLabel}${s}"
cp ../examples/Traj_${WindLabel}.ctr .
../bin/MetTraj_F Traj_${WindLabel}.ctr  > /dev/null 2>&1
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Error: MetTraj_F returned error code"
  exit 1
fi
###  Run check on ftraj1.dat and outputSonde2/ftraj1.dat
cmp -s ftraj1.dat output${WindLabel}2/ftraj1.dat
if [[ $? == 0 ]] ; then
  stat="PASS"
  printf " ---> ${GREEN}${stat}${NC}\n"
else
  stat="FAIL"
  printf " ---> ${RED}${stat}${NC}\n"
fi

sh clean.sh
echo "-----------------------------------------------------------"
echo "RUNNING SONDE TEST CASE 3: MetRegrid 2d"
echo "-----------------------------------------------------------"
s=1
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
cmp -s out02_1_2.dat output${WindLabel}3/out02_1_2.dat
if [[ $? == 0 ]] ; then
  stat="PASS"
  printf " ---> ${GREEN}${stat}${NC}\n"
else
  stat="FAIL"
  printf " ---> ${RED}${stat}${NC}\n"
fi


