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
    fileNames = cms.untracked.vstring(
      '/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00037C53-AAD1-E111-B1BE-003048D45F38.root', 
'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00050BBE-D5D2-E111-BB65-001E67398534.root', 
'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00B16DF1-8FD1-E111-ADBB-F04DA23BCE4C.root', 
'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00D370C2-E3D2-E111-9C87-003048673F0A.root', 
'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/00FB6563-58D2-E111-AFCB-001E67397D73.root', 
'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/023889B1-3ED2-E111-B9A4-00304866C51E.root', 
'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/026E25E5-65D2-E111-8481-00304867406E.root', 
'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/02A03572-32D2-E111-A90C-001E67397314.root', 
'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/02B0AB2A-ADD2-E111-8F7F-001E67397D7D.root', 
'/store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_S10_START53_V7A-v1/0000/02B1A617-E3D1-E111-B72F-001E67396E64.root'
    )
)

# Other statements
process.GlobalTag.globaltag = 'START52_V9::All'

process.Test = cms.EDAnalyzer("LHEPythiaEventAnalyzer",
    HistOutFile = cms.untracked.string('njetsMGOld.root'),
    theSrc = cms.untracked.string('source'),
    whichWeight = cms.untracked.int32(-1),
  #  isVBF = cms.untracked.bool(False),
  #  lookForGenProcID = cms.untracked.int32(0) # 24 = ZH, 26 = WH, 121,122 = ttH                         
)

process.p1 = cms.Path(process.Test)

