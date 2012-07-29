#! /bin/bash

if [ $# != 5 ] 
    then 
    echo oops...
    echo You are not using the right number of input arguments
    echo arg 1 - inputfile
    echo arg 2 - mZZlow
    echo arg 3 - mZZhigh
    echo arg 4 - true/false include PT
    echo arg 5 - true/false include Y
    echo Example: ./submitBatchJobs.sh inputfile 120 130 true true
    
    exit

    fi

infile=$1
mZZlow=$2
mZZhigh=$3
withPT=$4
withY=$5

root -l -n -q -b "Run_MELA.C(\"$infile\",${mZZlow},${mZZhigh},${withPT},${withY})"

echo "done"