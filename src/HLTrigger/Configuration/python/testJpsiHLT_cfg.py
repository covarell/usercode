import FWCore.ParameterSet.Config as cms

process = cms.Process("myHLT")

process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.VtxSmearedBetafuncEarlyCollision_cff")
process.load("Configuration.StandardSequences.Generator_cff")
process.load('Configuration/StandardSequences/GeometryPilot2_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load("Configuration.StandardSequences.L1Emulator_cff")
process.load("L1TriggerConfig.L1GtConfigProducers.Luminosity.lumi1030.L1Menu2008_2E30_cff")
process.load("HLTrigger.Configuration.HLT_cff")

process.schedule = process.HLTSchedule

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
     'rfio:/castor/cern.ch/user/c/covarell/Jpsi-nonp/inclJpsiee_GEN_SIM_RECO_11.root'
    )
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    makeTriggerResults = cms.untracked.bool(True)
)

process.out = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
    fileName = cms.untracked.string('file:jPsiHLTFromRaw.root'),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('HLT_DoubleEM8e29_Jpsi','HLT_DoubleEM8e29_Ups1s')
    )
)

process.pp = cms.EndPath(process.out)
process.schedule.append(process.pp)

process.GlobalTag.globaltag = 'STARTUP_V7::All'

