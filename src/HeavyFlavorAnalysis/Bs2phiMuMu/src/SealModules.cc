#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
// #include "HeavyFlavorAnalysis/Bs2phiMuMu/interface/Bs2phiMuMu.h"
#include "HeavyFlavorAnalysis/Bs2phiMuMu/interface/myAnalyzer.h"

DEFINE_SEAL_MODULE();
// DEFINE_ANOTHER_FWK_MODULE(Bs2phiMuMu);
DEFINE_ANOTHER_FWK_MODULE(myAnalyzer);
