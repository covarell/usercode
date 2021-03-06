import FWCore.ParameterSet.Config as cms

externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring('/cvmfs/cms.cern.ch/phys_generator/gridpacks/slc6_amd64_gcc700/13TeV/powheg/V2/ggHH_EWChL_NNPDF31_13TeV_M125_cHHH1/v3/ggHH_EWChL_slc6_amd64_gcc700_CMSSW_10_2_5_patch1_my_ggHH_EWChL.tgz'),
    nEvents = cms.untracked.uint32(5000),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
)

#Link to datacards:
#https://github.com/cms-sw/genproductions/tree/master/bin/MadGraph5_aMCatNLO/cards/production/2017/13TeV/exo_diboson/Spin-2/BulkGraviton_hh_narrow

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.MCTunes2017.PythiaCP5Settings_cfi import *
from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *
from Configuration.Generator.PSweightsPythia.PythiaPSweightsSettings_cfi import *


generator = cms.EDFilter("Pythia8HadronizerFilter",
                         maxEventsToPrint = cms.untracked.int32(1),
                         pythiaPylistVerbosity = cms.untracked.int32(1),
                         filterEfficiency = cms.untracked.double(1.0),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         comEnergy = cms.double(13000.),
                         PythiaParameters = cms.PSet(
                                                     pythia8CommonSettingsBlock,
                                                     pythia8CP5SettingsBlock,
                                                     pythia8PSweightsSettingsBlock,
                                                     pythia8PowhegEmissionVetoSettingsBlock,
                                                     processParameters = cms.vstring(
       'POWHEG:nFinal = 2',            
       '23:mMin = 0.05',
       '23:onMode = off',
       '23:onIfAny = 11 13 15', # only leptonic Z decays
       '25:m0 = 125.0',
       '25:onMode = off',
       '25:onIfMatch = 22 22',
       '25:onIfMatch = 23 23',
       'ResonanceDecayFilter:filter = on',
       'ResonanceDecayFilter:exclusive = on', #off: require at least the specified number of daughters, on: require exactly the specified number of daughters
       'ResonanceDecayFilter:eMuAsEquivalent = off', #on: treat electrons and muons as equivalent
       'ResonanceDecayFilter:eMuTauAsEquivalent = on', #on: treat electrons, muons , and taus as equivalent
       'ResonanceDecayFilter:allNuAsEquivalent = off', #on: treat all three neutrino flavours as equivalent
       'ResonanceDecayFilter:mothers = 25,23', #list of mothers not specified -> count all particles in hard process+resonance decays (better to avoid specifying mothers when including leptons from the lhe in counting, since intermediate resonances are not gauranteed to appear in general
       'ResonanceDecayFilter:daughters = 22,22,11,11,11,11',
                                                                                     ),
                                                     parameterSets = cms.vstring('pythia8CommonSettings',
                                                                                 'pythia8CP5Settings',
                                                                                 'pythia8PSweightsSettings',
 'pythia8PowhegEmissionVetoSettings',
                                                                                 'processParameters'
                                                                                 )
                                                     )
                         )


ProductionFilterSequence = cms.Sequence(generator)
