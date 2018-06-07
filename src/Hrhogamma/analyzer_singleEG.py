import FWCore.ParameterSet.Config as cms

process = cms.Process("REPORT")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 999999999
process.MessageLogger.categories.append('HLTrigReport')

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("RecoEgamma.PhotonIdentification.PhotonIDValueMapProducer_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

# trigger
import HLTrigger.HLTfilters.hltHighLevel_cfi

process.hltFilterMario = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterMario.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterMario.throw = cms.bool(False)
process.hltFilterMario.HLTPaths = ["HLT_Photon35_TwoProngs35_v*"]

process.triggerMario = cms.Path(process.hltFilterMario)

# PVs
process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)

# photons mediumID, no isolation yet
process.goodPhotons = cms.EDFilter("PATPhotonRefSelector",
   src = cms.InputTag("slimmedPhotons"),
   # cut = cms.string("pt>35 && (( abs(eta)<1.4442 && hadronicOverEm < 0.035 && sigmaIetaIeta < 0.0103) || ( 1.566<abs(eta)<2 && hadronicOverEm < 0.027 && sigmaIetaIeta < 0.0271))"),   ### Soffi
   cut = cms.string("pt>35 && r9 > 0.9 && hadronicOverEm < 0.05"),   
#   filter = cms.bool(True)
)

process.source = cms.Source("PoolSource",
     duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
      "file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_0.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_1.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_10.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_11.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_12.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_13.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_14.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_15.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_16.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_17.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_18.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_19.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_2.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_3.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_4.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_5.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_6.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_7.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_8.root",
"file:../../../../../Hrhogamma_singleEG/testMINIAOD_Hrhogamma_9.root",

        
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 10000 )
)

process.hlTrigReport = cms.EDAnalyzer("HLTrigReport",
    HLTriggerResults = cms.InputTag("TriggerResults","","HLT"),
    reportBy         = cms.untracked.string("job"),
    resetBy          = cms.untracked.string("never"),
    serviceBy        = cms.untracked.string("never")
)

process.analyzer = cms.EDAnalyzer("HMesonGammaAnalyzer",
    HistOutFile = cms.untracked.string('HrhogammaPlots_singleEG.root'),
    photonSrc = cms.InputTag('goodPhotons'),
    trackSrc = cms.InputTag('packedPFCandidates'),
    trackSrc2 = cms.InputTag('lostTracks'),
    HLTriggerResults = cms.InputTag("TriggerResults","","HLT"),
    HLTriggerObjects = cms.InputTag("slimmedPatTrigger"), 
    HLTriggerName = cms.string("HLT_Photon35_TwoProngs35_v1"), 
    hadronMass = cms.double(0.1396) #pi           0.4937  K              
)

process.mainPath = cms.Path(process.goodPrimaryVertices*process.goodPhotons*process.photonIDValueMapProducer*process.analyzer)
process.report = cms.EndPath( process.hlTrigReport )
