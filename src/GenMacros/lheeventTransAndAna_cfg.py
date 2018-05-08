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
    # 'file:/afs/cern.ch/work/c/covarell/powheg/POWHEG-BOX/VBF_H/testMSTW/pwgevents.lhe',    
#'/store/user/covarell//Gen/testMGOld/testMGOld_1.root',
#'/store/user/covarell//Gen/testMGOld/testMGOld_10.root',
#'/store/user/covarell//Gen/testMGOld/testMGOld_2.root',
#'/store/user/covarell//Gen/testMGOld/testMGOld_4.root',
#'/store/user/covarell//Gen/testMGOld/testMGOld_5.root',
#'/store/user/covarell//Gen/testMGOld/testMGOld_6.root',
#'/store/user/covarell//Gen/testMGOld/testMGOld_7.root',
#'/store/user/covarell//Gen/testMGOld/testMGOld_8.root',
#'/store/user/covarell//Gen/testMGOld/testMGOld_9.root',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10001.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10002.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10003.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10004.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10005.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10006.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10007.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10008.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10009.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10010.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10011.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10012.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10013.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10014.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10015.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10016.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10017.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10018.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10019.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10020.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10021.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10022.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10023.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10024.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10025.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10026.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10027.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10028.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10029.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10030.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10031.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10032.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10033.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10034.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10035.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10036.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10037.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10038.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10039.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10040.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10041.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10042.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10043.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10044.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10045.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10046.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10047.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10048.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10049.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10050.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10051.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10052.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10053.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10054.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10055.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10056.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10057.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10058.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10059.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10060.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10061.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10062.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10063.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10064.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10065.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10066.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10067.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10068.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10069.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10070.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10071.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10072.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10073.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10074.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10075.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10076.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10077.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10078.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10079.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10080.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10081.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10082.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10083.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10084.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10085.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10086.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10087.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10088.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10089.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10090.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10091.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10092.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10093.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10094.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10095.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10096.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10097.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10098.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10099.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10100.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10101.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10102.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10103.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10104.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10105.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10106.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10107.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10108.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10109.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10110.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10111.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10112.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10113.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10114.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10115.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10116.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10117.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10118.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10119.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10120.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10121.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10123.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10124.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10125.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10126.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10127.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10128.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10129.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10130.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10131.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10132.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10133.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10134.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10135.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10136.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10137.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10138.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10139.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10140.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10141.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10142.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10143.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10144.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10145.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10146.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10147.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10148.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10149.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10150.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10151.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10152.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10153.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10154.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10155.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10156.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10157.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10158.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10159.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10160.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10161.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10162.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10163.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10164.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10165.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10166.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10167.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10168.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10169.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10170.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10171.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10172.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10173.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10174.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10175.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10176.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10177.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10178.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10179.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10180.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10181.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10182.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10183.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10184.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10185.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10186.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10187.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10188.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10189.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10190.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10191.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10192.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10193.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10194.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10195.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10196.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10197.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10198.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10199.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10200.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10201.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10202.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10203.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10204.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10205.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10206.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10207.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10208.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10209.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10210.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10211.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10212.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10213.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10214.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10215.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10216.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10217.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10218.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10219.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10220.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10221.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10222.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10223.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10224.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10225.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10226.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10227.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10228.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10229.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10230.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10231.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10232.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10233.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10234.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10235.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10236.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10237.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10238.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10239.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10240.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10241.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10242.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10243.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10244.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10245.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10246.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10247.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10248.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10249.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10250.lhe',
         '/store/lhe/5591/DYJetsToLL_M-50_8TeV-madgraph-tarball_10251.lhe'

   )
)

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.Test = cms.EDAnalyzer("LHEEventAnalyzer",
    HistOutFile = cms.untracked.string('multMGOld.root'),
    theSrc = cms.untracked.string('source')                           
)

process.p1 = cms.Path(process.Test)

