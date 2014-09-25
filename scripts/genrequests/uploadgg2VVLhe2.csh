#!/usr/bin/tcsh

# Change process name in Matthias' config files
rm *.lhe lista*
~/myscripts/cmsStageEOSDir.pl /store/user/covarell/GG2VV-lhe/${1} lhe tempor1 200
~/myscripts/cmsStageEOSDir.pl /store/user/covarell/GG2VV-lhe/${1}2 lhe tempor2 200
~/myscripts/cmsStageEOSDir.pl /store/user/covarell/GG2VV-lhe/${1}3 lhe tempor3 200
~/myscripts/cmsStageEOSDir.pl /store/user/covarell/GG2VV-lhe/${1}4 lhe tempor2 200
~/myscripts/cmsStageEOSDir.pl /store/user/covarell/GG2VV-lhe/${1}5 lhe tempor3 200
#~/myscripts/cmsStageEOSDir.pl /store/user/covarell/GG2VV-lhe/${1}6 lhe tempor2 200
#~/myscripts/cmsStageEOSDir.pl /store/user/covarell/GG2VV-lhe/${1}7 lhe tempor3 200
#~/myscripts/cmsStageEOSDir.pl /store/user/covarell/GG2VV-lhe/${1}8 lhe tempor2 200

foreach nummero (0 1 2 3 4 5 6 7 8 9)
  ls ./*${nummero}.lhe > lista${nummero}
  ~/work/gg2VV/gg2VV-3.1.5/gg2VV/modifyLheHeadings lista${nummero}
  mv out.lhe $1_tot${nummero}.lhe
end

/afs/cern.ch/user/c/covarell/scratch0/mcrequests/CMSSW_5_3_3/src/GeneratorInterface/LHEInterface/scripts/cmsLHEtoEOSManager.py -n 8000 -f $1_tot0.lhe,$1_tot1.lhe,$1_tot2.lhe,$1_tot3.lhe,$1_tot4.lhe,$1_tot5.lhe,$1_tot6.lhe,$1_tot7.lhe,$1_tot8.lhe,$1_tot9.lhe

#echo Process name changed from $1 to $2 in all files
 
