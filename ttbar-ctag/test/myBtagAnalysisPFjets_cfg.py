# The following comments couldn't be translated into the new config version:

#! /bin/env cmsRun

import FWCore.ParameterSet.Config as cms

process = cms.Process("analysis")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START53_V23::All' 
process.load("DQMServices.Components.DQMEnvironment_cfi")

#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("DQMServices.Core.DQM_cfg")

process.load("RecoBTag.Configuration.RecoBTag_cff")

process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")

process.load("Validation.RecoB.rerunBtag_cfi")  

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso
process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')

from DQMOffline.RecoB.bTagCommon_cff import *
from DQMOffline.RecoB.bTagCombinedSVAnalysis_cff import *

process.bTagValidation = cms.EDAnalyzer("myBTagAnalyzerMC",
    bTagCommonBlock,
    jetCorrection = cms.string(''),
    recJetMatching = cms.PSet(
        refJetCorrection = cms.string(''),
        recJetCorrection = cms.string(''),
        maxChi2 = cms.double(50),
        # Corrected calo jets
        sigmaDeltaR = cms.double(0.1),
        sigmaDeltaE = cms.double(0.15)
    ),
    tagConfig = cms.PSet(
        bTagCombinedSVAnalysisBlock,
        tagInfoLabel1 = cms.InputTag("newImpactParameterTagInfos"),
        tagInfoLabel2 = cms.InputTag("newSecondaryVertexTagInfos"),
        tagInfoLabel3 = cms.InputTag("none"),  
        tagLabel = cms.InputTag("newCombinedSecondaryVertexBJetTags")
    ), 
    vtxColl = cms.InputTag("offlinePrimaryVertices"),
    electronColl = cms.InputTag("gsfElectrons"),
    muonColl = cms.InputTag("muons"),
    metColl = cms.InputTag("pfMet"),
    jetColl = cms.InputTag("patJets"),                                    
    NeventsTOT = cms.int32(9124),
    xsec = cms.double(225.197),
    lumi = cms.double(19702.),
)
process.bTagValidation.jetMCSrc = 'AK5byValAlgo'

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

process.plots = cms.Path( process.pfParticleSelectionSequence *
                          process.eleIsoSequence *
                          process.newJetTracksAssociator *
                          process.newJetBtagging *
                          process.myPartons *
                          process.AK5Flavour *
                          process.bTagValidation)
#process.dqmEnv.subSystemFolder = 'BTAG'
#process.dqmSaver.producer = 'DQM'
#process.dqmSaver.workflow = '/POG/BTAG/BJET'
#process.dqmSaver.convention = 'Offline'
#process.dqmSaver.saveByRun = cms.untracked.int32(-1)
#process.dqmSaver.saveAtJobEnd =cms.untracked.bool(True) 
#process.dqmSaver.forceRunNumber = cms.untracked.int32(1)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('ttbar_ntuple_test.root')
)

process.PoolSource.fileNames = [
'/store/relval/CMSSW_5_3_14/RelValTTbar/GEN-SIM-RECO/START53_LV4_Feb7-v2/00000/7E099684-8390-E311-96F7-003048FFCC0A.root',
'/store/relval/CMSSW_5_3_14/RelValTTbar/GEN-SIM-RECO/START53_LV4_Feb7-v2/00000/F89412D1-8790-E311-8799-0025905A6070.root'
]

