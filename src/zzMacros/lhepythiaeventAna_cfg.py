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
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    # 'file:/afs/cern.ch/work/c/covarell/powheg/POWHEG-BOX/VBF_H/testMSTW/pwgevents.lhe',    
    'file:/tmp/covarell/POWHEG_PYTHIA6_Tauola_H125_ZZ_4l_8TeV_noUEnoHadr_cff_py_GEN.root',
    
   
    )
)

# Other statements
process.GlobalTag.globaltag = 'START52_V9::All'

process.Test = cms.EDAnalyzer("LHEPythiaEventAnalyzer",
    HistOutFile = cms.untracked.string('gg400_pythia.root'),
    theSrc = cms.untracked.string('source'),
    isVBF = cms.untracked.bool(False),
    lookForGenProcID = cms.untracked.int32(0) # 24 = ZH, 26 = WH, 121,122 = ttH                         
)

process.p1 = cms.Path(process.Test)

