#! /bin/bash

if [ $# != 6 ] 
    then 
    echo oops...
    echo You are not using the right number of input arguments
    echo arg 1 - inputfile
    echo arg 2 - energy 7 or 8 TeV
    echo arg 3 - mZZlow
    echo arg 4 - mZZhigh
    echo arg 5 - true/false include PT
    echo arg 6 - true/false include Y
    echo Example: ./submitBatchJobs.sh inputfile 120 130 true true
    
    exit

    fi

infile=$1
energy=$2
mZZlow=$3
mZZhigh=$4
withPT=$5
withY=$6

root -l -n -q -b "Run_MELA.C(\"$infile\",${energy},${mZZlow},${mZZhigh},${withPT},${withY})"

echo "done"