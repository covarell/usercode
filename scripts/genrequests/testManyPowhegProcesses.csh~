#!/usr/bin/tcsh

# Change process name in Matthias' config files
rm *.lhe lista*

foreach nummero (0 1 2 3 4 5 6 7 8 9)
  ls ~mdalchen/public/LHEfiles/$1/gg2VV-*${nummero}_unweightedEvents.dat > lista${nummero}
  ~/work/gg2VV/gg2VV-3.1.5/gg2VV/modifyLheHeadings lista${nummero}
  mv out.lhe $1_${nummero}.lhe
end

/afs/cern.ch/user/c/covarell/scratch0/mcrequests/CMSSW_5_3_3/src/GeneratorInterface/LHEInterface/scripts/cmsLHEtoEOSManager.py -n 8000 -f $1_0.lhe,$1_1.lhe,$1_2.lhe,$1_3.lhe,$1_4.lhe,$1_5.lhe,$1_6.lhe,$1_7.lhe,$1_8.lhe,$1_9.lhe

#echo Process name changed from $1 to $2 in all files
 
