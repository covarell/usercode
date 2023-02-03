import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500000)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# Input source
process.source = cms.Source("LHESource",
  #  duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/c/covarell/mcfm_klambda/ggZZ_2e2mu_c6eq0_MCFM8/cmsgrid_final.lhe',    

   )
)

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.Test = cms.EDAnalyzer("LHEEventAnalyzer",
    HistOutFile = cms.untracked.string('c6eq0.root'),
    theSrc = cms.untracked.string('source'),                           
    #xSection = cms.untracked.double(2.0807),
    xSection = cms.untracked.double(1.7265),
)

process.p1 = cms.Path(process.Test)

