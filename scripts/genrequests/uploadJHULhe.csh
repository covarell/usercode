#!/usr/bin/tcsh

# Change process name in Matthias' config files

foreach nummero (140 145 175)
  rm *.lhe lista    
  cp /tmp/joaguila/$nummero/*.lhe .
  ls *.lhe > lista
  mergeLheFiles lista
  mv out.lhe SMHiggsToZZTo4L_M-${nummero}_7TeV_POWHEG-JHUgenV3-pythia6.lhe
  /afs/cern.ch/user/c/covarell/scratch0/mcrequests/CMSSW_5_3_3/src/GeneratorInterface/LHEInterface/scripts/cmsLHEtoEOSManager.py -n 8000 -f SMHiggsToZZTo4L_M-${nummero}_7TeV_POWHEG-JHUgenV3-pythia6.lhe
end

#echo Process name changed from $1 to $2 in all files
 
