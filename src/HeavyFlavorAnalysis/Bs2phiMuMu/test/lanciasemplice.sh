#!/bin/sh

source /opt/exp_soft/cms/cmsset_default.sh

cd /data06/users/covarell/CMSSW_1_6_11/src/HeavyFlavorAnalysis/Bs2phiMuMu/test
# eval \`scramv1 runtime -sh\`
source /data06/users/covarell/CMSSW_1_6_11/src/HeavyFlavorAnalysis/Bs2phiMuMu/test/patchScramRuntime.sh

cmsRun BsphiMuMu_SIM_DIGI_RECO_noinput.cfg >& /data06/users/covarell/cmssw/BsphiMuMu_SIM_DIGI_RECO.log
touch DONE

