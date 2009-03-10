#!/bin/sh

export MYIFILE=$1
export MYOFILE=$2

cat << eof >& tmpcfg/tmp.cfg

source = PoolSource { 		
   	untracked vstring fileNames = { 
         "file:/home/covarell/home-data06/cmssw/$MYIFILE.root"
        }	
    }

module FEVT = PoolOutputModule { 
    using FEVTSIMEventContent
    untracked string fileName = "/home/covarell/home-data06/cmssw/tmpntpl/$MYOFILE.root"
    untracked PSet dataset = {	
      untracked string dataTier = "GEN-SIM-DIGI-RECO"
    }   
   }

}
eof

cp BsphiMuMu_SIM_DIGI_RECO_noinput.cfg tmpcfg/$MYIFILE.cfg
cat tmpcfg/tmp.cfg >> tmpcfg/$MYIFILE.cfg
# rm -f prova.cfg

# Generate submitting script
cat << eof >& tmpscripts/$MYIFILE.sh
source /opt/exp_soft/cms/cmsset_default.sh

cd /data06/users/covarell/CMSSW_1_6_11/src/HeavyFlavorAnalysis/Bs2phiMuMu/test
# eval \`scramv1 runtime -sh\`
source /data06/users/covarell/CMSSW_1_6_11/src/HeavyFlavorAnalysis/Bs2phiMuMu/test/patchScramRuntime.sh

cmsRun tmpcfg/$MYIFILE.cfg >& /data06/users/covarell/cmssw/tmpntpl/$MYOFILE.log
touch DONE
eof
 
# Submit job
chmod ugo+x tmpscripts/$MYIFILE.sh
qsub -q babar-infinite tmpscripts/$MYIFILE.sh


