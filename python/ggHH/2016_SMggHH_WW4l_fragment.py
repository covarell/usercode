import FWCore.ParameterSet.Config as cms

# link to cards:
# https://github.com/cms-sw/genproductions/tree/dac0b0a9a49aeaddc656529d275e5d426effce44/bin/MadGraph5_aMCatNLO/cards/production/13TeV/NonRes_hh/GF_HH_node_SM
# decay fragment from example:
# https://github.com/cms-sw/genproductions/blob/2e8edbb4c940bf3eea6c3f7af51727a6eb545d8a/python/ThirteenTeV/Higgs/HH/ResonanceDecayFilter_example_HHTo2B2L2Nu_madgraph_pythia8_cff.py


externalLHEProducer = cms.EDProducer("ExternalLHEProducer",
    args = cms.vstring('/cvmfs/cms.cern.ch/phys_generator/gridpacks/slc6_amd64_gcc700/13TeV/powheg/V2/ggHH_EWChL_NNPDF31_13TeV_M125_cHHH1/v3/ggHH_EWChL_slc6_amd64_gcc700_CMSSW_10_2_5_patch1_my_ggHH_EWChL.tgz'),
    nEvents = cms.untracked.uint32(5000),
    numberOfParameters = cms.uint32(1),
    outputFile = cms.string('cmsgrid_final.lhe'),
    scriptName = cms.FileInPath('GeneratorInterface/LHEInterface/data/run_generic_tarball_cvmfs.sh')
)

from Configuration.Generator.Pythia8CommonSettings_cfi import *
from Configuration.Generator.Pythia8CUEP8M1Settings_cfi import *
from Configuration.Generator.Pythia8PowhegEmissionVetoSettings_cfi import *

generator = cms.EDFilter("Pythia8HadronizerFilter",
                         maxEventsToPrint = cms.untracked.int32(1),
                         pythiaPylistVerbosity = cms.untracked.int32(1),
                         filterEfficiency = cms.untracked.double(1.0),
                         pythiaHepMCVerbosity = cms.untracked.bool(False),
                         comEnergy = cms.double(13000.),
                         PythiaParameters = cms.PSet(
        pythia8CommonSettingsBlock,
        pythia8CUEP8M1SettingsBlock,
        pythia8PowhegEmissionVetoSettingsBlock,
        processParameters = cms.vstring(
            'POWHEG:nFinal = 2', 
            '23:mMin = 0.05',
            '23:onMode = off',
            '23:onIfAny = 11 13 15', # only leptonic Z decays
            '25:m0 = 125.0',
            '25:onMode = off',
            '25:onIfMatch = 24 24',
            '25:onIfMatch = 23 23',
            'ResonanceDecayFilter:filter = on',
            'ResonanceDecayFilter:exclusive = off', #off: require at least the specified number of daughters, on: require exactly the specified number of daughters
            'ResonanceDecayFilter:eMuAsEquivalent = off', #on: treat electrons and muons as equivalent
            'ResonanceDecayFilter:eMuTauAsEquivalent = on', #on: treat electrons, muons , and taus as equivalent
            'ResonanceDecayFilter:allNuAsEquivalent = off', #on: treat all three neutrino flavours as equivalent
            'ResonanceDecayFilter:mothers = 23,24', # list of mothers not specified -> count all particles in hard process+resonance decays (better to avoid specifying mothers when including leptons from the lhe in counting, since intermediate resonances are not gauranteed to appear in general
            'ResonanceDecayFilter:daughters = 11,11,11,11,11',
            ),
        parameterSets = cms.vstring('pythia8CommonSettings',
                                    'pythia8CUEP8M1Settings',
                                    'pythia8PowhegEmissionVetoSettings', 
                                    'processParameters'
                                    )
        )
                         )


ProductionFilterSequence = cms.Sequence(generator)
