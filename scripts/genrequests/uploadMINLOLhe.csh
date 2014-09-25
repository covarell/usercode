#!/usr/bin/tcsh

# Change process name in Matthias' config files

rm * 
/afs/cern.ch/user/c/covarell/myscripts/cmsStageEOSDir.pl /store/cmst3/user/govoni/powheg/prodU04_$1 lhe boh 50
/afs/cern.ch/user/c/covarell/myscripts/genRequests/removePOWHEGTags.pl

/afs/cern.ch/user/c/covarell/scratch0/mcrequests/CMSSW_5_3_3/src/GeneratorInterface/LHEInterface/scripts/cmsLHEtoEOSManager.py -n 8000 -f pwgevents-0001.lhe,pwgevents-0002.lhe,pwgevents-0003.lhe,pwgevents-0004.lhe,pwgevents-0005.lhe,pwgevents-0006.lhe,pwgevents-0007.lhe,pwgevents-0008.lhe,pwgevents-0009.lhe,pwgevents-0010.lhe > tempofile 

grep 'Creating' tempofile > bookkeep$1.list
rm tempofile

/afs/cern.ch/user/c/covarell/scratch0/mcrequests/CMSSW_5_3_3/src/GeneratorInterface/LHEInterface/scripts/cmsLHEtoEOSManager.py -n 8000 -f pwgevents-0011.lhe,pwgevents-0012.lhe,pwgevents-0013.lhe,pwgevents-0014.lhe,pwgevents-0015.lhe,pwgevents-0016.lhe,pwgevents-0017.lhe,pwgevents-0018.lhe,pwgevents-0019.lhe,pwgevents-0020.lhe > tempofile 

grep 'Creating' tempofile > bookkeep$1_2.list
rm tempofile

/afs/cern.ch/user/c/covarell/scratch0/mcrequests/CMSSW_5_3_3/src/GeneratorInterface/LHEInterface/scripts/cmsLHEtoEOSManager.py -n 8000 -f pwgevents-0021.lhe,pwgevents-0022.lhe,pwgevents-0023.lhe,pwgevents-0024.lhe,pwgevents-0025.lhe,pwgevents-0026.lhe,pwgevents-0027.lhe,pwgevents-0028.lhe,pwgevents-0029.lhe,pwgevents-0030.lhe > tempofile 

grep 'Creating' tempofile > bookkeep$1_3.list
rm tempofile

/afs/cern.ch/user/c/covarell/scratch0/mcrequests/CMSSW_5_3_3/src/GeneratorInterface/LHEInterface/scripts/cmsLHEtoEOSManager.py -n 8000 -f pwgevents-0031.lhe,pwgevents-0032.lhe,pwgevents-0033.lhe,pwgevents-0034.lhe,pwgevents-0035.lhe,pwgevents-0036.lhe,pwgevents-0037.lhe,pwgevents-0038.lhe,pwgevents-0039.lhe,pwgevents-0040.lhe > tempofile 

grep 'Creating' tempofile > bookkeep$1_4.list
rm tempofile

/afs/cern.ch/user/c/covarell/scratch0/mcrequests/CMSSW_5_3_3/src/GeneratorInterface/LHEInterface/scripts/cmsLHEtoEOSManager.py -n 8000 -f pwgevents-0041.lhe,pwgevents-0042.lhe,pwgevents-0043.lhe,pwgevents-0044.lhe,pwgevents-0045.lhe,pwgevents-0046.lhe,pwgevents-0047.lhe,pwgevents-0048.lhe,pwgevents-0049.lhe,pwgevents-0050.lhe > tempofile 

grep 'Creating' tempofile > bookkeep$1_5.list
rm tempofile

/afs/cern.ch/user/c/covarell/myscripts/genRequests/bookkeepMINLO.pl $1

#echo Process name changed from $1 to $2 in all files
 
