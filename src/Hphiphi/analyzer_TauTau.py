import FWCore.ParameterSet.Config as cms

process = cms.Process("REPORT")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.categories.append('HLTrigReport')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
#process.load("RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

# trigger
#import HLTrigger.HLTfilters.hltHighLevel_cfi

#process.hltFilterMario = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
#process.hltFilterMario.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
#process.hltFilterMario.throw = cms.bool(False)
#process.hltFilterMario.HLTPaths = ["HLT_Photon35_TwoProngs35_v*"]

#process.triggerMario = cms.Path(process.hltFilterMario)

# PVs
process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)

# photons mediumID, no isolation yet
#process.goodPhotons = cms.EDFilter("PATPhotonRefSelector",
#   src = cms.InputTag("slimmedPhotons"),
 #  cut = cms.string("pt>35 && (( abs(eta)<1.4442 && hadronicOverEm < 0.035 && sigmaIetaIeta < 0.0103) || ( 1.566<abs(eta)<2 && hadronicOverEm < 0.027 && sigmaIetaIeta < 0.0271))"),   ### Soffi
#   cut = cms.string("pt>35 && abs(eta) < 2.1 && r9 > 0.9 && hadronicOverEm < 0.05"),  
#   filter = cms.bool(True)
#)

process.source = cms.Source("PoolSource",
     duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
        "file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_0.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_1.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_10.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_11.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_12.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_13.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_14.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_15.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_16.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_17.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_18.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_19.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_2.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_20.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_21.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_22.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_23.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_24.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_25.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_26.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_27.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_28.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_29.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_3.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_30.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_31.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_32.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_33.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_34.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_35.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_36.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_37.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_38.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_39.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_4.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_40.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_41.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_42.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_43.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_44.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_45.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_46.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_47.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_48.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_49.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_5.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_6.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_7.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_8.root",
"file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphiphi/HphiphiEvtGen_9.root",

    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 100000 )
)

process.hlTrigReport = cms.EDAnalyzer("HLTrigReport",
    HLTriggerResults = cms.InputTag("TriggerResults","","HLT"),
    reportBy         = cms.untracked.string("job"),
    resetBy          = cms.untracked.string("never"),
    serviceBy        = cms.untracked.string("never")
)

process.analyzer = cms.EDAnalyzer("HPhiPhiAnalyzer",
    HistOutFile = cms.untracked.string('HphiphiPlots.root'),
    #photonSrc = cms.InputTag('goodPhotons'),
    trackSrc = cms.InputTag('packedPFCandidates'),
    trackSrc2 = cms.InputTag('lostTracks'),
    HLTriggerResults = cms.InputTag("TriggerResults","","HLT"),
    HLTriggerObjects = cms.InputTag("slimmedPatTrigger"),     
    HLTriggerName = cms.string("HLT_DoublePNetTauhPFJet30_Medium_L2NN_eta2p3_v*"), 
    hadronMass = cms.double(0.4937) #pi           0.4937  K              
)

process.mainPath = cms.Path(process.goodPrimaryVertices*process.goodPhotons*process.photonIDValueMapProducer*process.analyzer)
process.report = cms.EndPath( process.hlTrigReport )
