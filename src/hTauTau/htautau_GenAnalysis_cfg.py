import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")
process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        )
    ),
    destinations = cms.untracked.vstring('cout')
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/c/covarell/cmst3/hTauTauTest.root')
    # fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/c/covarell/cmst3/SUSYhTauTauTest.root')
    # fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/c/covarell/cmst3/zTauTauTest.root')
)

process.Test = cms.EDAnalyzer("HTauTauGenAnalyzer",
    HistOutFile = cms.untracked.string('testSMHiggs.root'),
    # HistOutFile = cms.untracked.string('testSUSYHiggs.root'),
    # HistOutFile = cms.untracked.string('testZ.root'),
    theSrc = cms.untracked.string('generator')
)

process.p1 = cms.Path(process.Test)

