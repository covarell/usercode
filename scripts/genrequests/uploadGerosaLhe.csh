#!/usr/bin/tcsh

# Change process name in Matthias' config files

cmsStage -f /store/caf/user/rgerosa/powheg_lhe_production/ggH_CPS/$1.tar.gz .
gunzip $1.tar.gz
tar -xvf $1.tar
cp pwgevents.lhe gg_8TeV_$2.lhe
/afs/cern.ch/user/c/covarell/myscripts/genRequests/removePOWHEGTags.pl
/afs/cern.ch/user/c/covarell/scratch0/mcrequests/CMSSW_5_3_3/src/GeneratorInterface/LHEInterface/scripts/cmsLHEtoEOSManager.py -n 8000 -f gg_8TeV_$2.lhe
rm -f gg_8TeV_$2.lhe joboutput.txt

echo Process name changed from $1 to $2 in all files
 
