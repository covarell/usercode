# Auto generated configuration file
# using: 
# Revision: 1.381.2.6 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: MCDBtoEDM --conditions START50_V12::All -s NONE --eventcontent RAWSIM --datatier GEN --filein file:/tmp/covarell/h_PG_TT_TTVBF_90_0.lhe --no_exec -n 1
import FWCore.ParameterSet.Config as cms

process = cms.Process('LHE')

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("LHESource",
    fileNames = cms.untracked.vstring(
                                      'file:/tmp/covarell/h_PG_TT_TTVBFnew_90_0.lhe',
                                      'file:/tmp/covarell/h_PG_TT_TTVBFnew_90_1.lhe',
                                      'file:/tmp/covarell/h_PG_TT_TTVBFnew_90_2.lhe',
                                      'file:/tmp/covarell/h_PG_TT_TTVBFnew_90_3.lhe',
                                      'file:/tmp/covarell/h_PG_TT_TTVBFnew_90_4.lhe',
                                      # 'file:/tmp/covarell/h_PG_TT_TT_90_0.lhe'   
                                      )
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.381.2.6 $'),
    annotation = cms.untracked.string('MCDBtoEDM nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('MCDBtoEDM_NONE.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN')
    )
)

# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START50_V12::All'

# Path and EndPath definitions
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.RAWSIMoutput_step)

def customise(process):

    #Tweak Message logger to dump G4cout and G4cerr messages in G4msg.log
    #print process.MessageLogger.__dict__
    process.MessageLogger.destinations=cms.untracked.vstring('cout'
                                                             ,'cerr'
                                                             )
    process.MessageLogger.categories=cms.untracked.vstring('FwkJob'
                                                           ,'FwkReport'
                                                           ,'FwkSummary'
                                                           ,'Root_NoDictionary'
                                                           ,'Generator'
                                                           ,'LHEInterface'
                                                           )
    #Configuring the Gen.log output
    process.MessageLogger.cerr =  cms.untracked.PSet(
        noTimeStamps = cms.untracked.bool(True)
        #First eliminate unneeded output
        ,threshold = cms.untracked.string('INFO')
        ,INFO = cms.untracked.PSet(limit = cms.untracked.int32(-1))
        ,FwkReport = cms.untracked.PSet(limit = cms.untracked.int32(-1),
                                   reportEvery = cms.untracked.int32(1000))
        ,FwkSummary = cms.untracked.PSet(limit = cms.untracked.int32(-1))
        ,Root_NoDictionary = cms.untracked.PSet(limit = cms.untracked.int32(-1))
        ,FwkJob = cms.untracked.PSet(limit = cms.untracked.int32(-1))
        ,Generator = cms.untracked.PSet(limit = cms.untracked.int32(0))
        ,LHEInterface = cms.untracked.PSet(limit = cms.untracked.int32(10000))
        )

    #Add these 3 lines to put back the summary for timing information at the end of the logfile
    #(needed for TimeReport report)
    process.options = cms.untracked.PSet(
        wantSummary = cms.untracked.bool(True)
        )

        
    return(process)

process = customise(process) 
