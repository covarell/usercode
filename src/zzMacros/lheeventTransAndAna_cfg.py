import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger = cms.Service("MessageLogger",
    cout = cms.untracked.PSet(
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        )
    ),
    destinations = cms.untracked.vstring('cout')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring(
    # 'file:/afs/cern.ch/work/c/covarell/powheg/POWHEG-BOX/VBF_H/testMSTW/pwgevents.lhe',    
    'file:/afs/cern.ch/work/c/covarell/powheg/POWHEG-BOX/gg_H_quark-mass-effects/testrun-lhc/pwgevents200.lhe',
   
    )
)

# Other statements
process.GlobalTag.globaltag = 'START50_V13::All'

process.Test = cms.EDAnalyzer("LHEEventAnalyzer",
    HistOutFile = cms.untracked.string('ggH200_finiteMT.root'),
    theSrc = cms.untracked.string('source'),
    isVBF = cms.untracked.bool(False)                           
)

process.p1 = cms.Path(process.Test)

