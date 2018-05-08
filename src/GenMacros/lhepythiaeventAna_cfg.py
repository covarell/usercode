import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYSIS")

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(500000)
)

# Input source
process.source = cms.Source("PoolSource",
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_100.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_101.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_102.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_103.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_104.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_105.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_106.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_107.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_108.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_109.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_110.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_111.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_112.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_113.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_114.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_115.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_116.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_117.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_118.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_119.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_120.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_121.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_122.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_123.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_124.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_125.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_126.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_127.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_128.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_129.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_130.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_131.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_132.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_133.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_134.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_135.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_136.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_137.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_138.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_139.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_140.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_141.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_142.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_143.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_144.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_145.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_146.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_147.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_148.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_149.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_150.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_151.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_152.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_153.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_154.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_155.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_156.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_157.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_158.root',
     'file:../../../../../HWJvalid/testGEN_HWminusJ_WMuNu_159.root',

    )
)

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.Test = cms.EDAnalyzer("LHEPythiaEventAnalyzer",
    HistOutFile = cms.untracked.string('HWJplots.root'),
    theSrc = cms.untracked.string('externalLHEProducer'),
    whichWeight = cms.untracked.int32(-1),
  #  isVBF = cms.untracked.bool(False),
  #  lookForGenProcID = cms.untracked.int32(0) # 24 = ZH, 26 = WH, 121,122 = ttH                         
)

process.p1 = cms.Path(process.Test)

