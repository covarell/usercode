#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "GeneratorInterface/EvtGenInterface/interface/EvtGenProducer.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(EvtGenProducer);

