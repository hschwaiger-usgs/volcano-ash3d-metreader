#!/bin/bash

# This bash script is run for the target 'make check'

rc=0
if [ -z ${WINDROOT} ];then
 # Standard Linux location
 WINDROOT="/data/WindFiles"
 # Mac
 #WINDROOT="/opt/data/WindFiles"
fi

# Check to see if the NCEP data for 1980 is present
ls -1r ${WINDROOT}/NCEP/1980/air.1980.nc
rc=$((rc + $?))
if [[ "$rc" -gt 0 ]] ; then
  echo "Warning: Could not find NCEP data for 1980"
  echo "Only the Sonde test cases will be run."
  echo "To download the NCEP data, run:"
  echo "/opt/USGS/bin/autorun_scripts/autorun_scripts/get_NCEP_50YearReanalysis.sh 1980"
  RUNNCEP=F
  nwin=1
else
  RUNNCEP=T
  nwin=2
fi

pushd tests
bash clean.sh
bash run_tests_Sonde.sh
bash clean.sh
if [[ "$RUNNCEP" -eq T ]] ; then
  bash run_tests_NCEP.sh
  bash clean.sh
fi
popd
