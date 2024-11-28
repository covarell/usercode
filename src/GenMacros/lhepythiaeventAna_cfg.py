import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500000)
)

# Input source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_0.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_10.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_11.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_13.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_14.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_15.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_17.root',
#'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_19.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_1.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_2.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_3.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_4.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_5.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_6.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_8.root',
'file:/eos/user/p/pellicci/MesonGamma_root/2024/testOfficialProd/ggHphigamma/HphigammaEvtGen_9.root',

    )
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.Test = cms.EDAnalyzer("LHEPythiaEventAnalyzer",
    HistOutFile = cms.untracked.string('HphigHelPlots.root'),
    theSrc = cms.untracked.string('externalLHEProducer'),
    whichWeight = cms.untracked.int32(-1),
  #  isVBF = cms.untracked.bool(False),
  #  lookForGenProcID = cms.untracked.int32(0) # 24 = ZH, 26 = WH, 121,122 = ttH                         
)

process.p1 = cms.Path(process.Test)

