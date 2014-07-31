import FWCore.ParameterSet.Config as cms

from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import *
from RecoBTag.Configuration.RecoBTag_cff import *

# create a new jets and tracks association
newJetTracksAssociatorAtVertex = ic5JetTracksAssociatorAtVertex.clone()
newJetTracksAssociatorAtVertex.jets = "ak5PFJets"
newJetTracksAssociatorAtVertex.tracks = "generalTracks"

# impact parameter b-tag
newImpactParameterTagInfos = impactParameterTagInfos.clone()
newImpactParameterTagInfos.jetTracks = "newJetTracksAssociatorAtVertex"
newTrackCountingHighEffBJetTags = trackCountingHighEffBJetTags.clone()
newTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
newTrackCountingHighPurBJetTags = trackCountingHighPurBJetTags.clone()
newTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
newJetProbabilityBJetTags = jetProbabilityBJetTags.clone()
newJetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )
newJetBProbabilityBJetTags = jetBProbabilityBJetTags.clone()
newJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos") )

newSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
newSecondaryVertexTagInfos.trackIPTagInfos = "newImpactParameterTagInfos"
newSimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags.clone()
newSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSecondaryVertexTagInfos") )
newCombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
newCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"), cms.InputTag("newSecondaryVertexTagInfos") )
newCombinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTags.clone()
newCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("newImpactParameterTagInfos"), cms.InputTag("newSecondaryVertexTagInfos") )

# soft electron b-tag
newSoftPFElectronsTagInfos = softPFElectronsTagInfos.clone()
newSoftPFElectronsTagInfos.jets = "ak5PFJets"
newSoftPFElectronBJetTags = softPFElectronBJetTags.clone()
newSoftPFElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftPFElectronsTagInfos") )

# soft muon b-tag
newSoftPFMuonsTagInfos = softPFMuonsTagInfos.clone()
newSoftPFMuonsTagInfos.jets = "ak5PFJets"
newSoftPFMuonBJetTags = softPFMuonBJetTags.clone()
newSoftPFMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftPFMuonsTagInfos") )
newSoftPFMuonByIP3dBJetTags = softPFMuonByIP3dBJetTags.clone()
newSoftPFMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftPFMuonsTagInfos") )
newSoftPFMuonByPtBJetTags = softPFMuonByPtBJetTags.clone()
newSoftPFMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("newSoftPFMuonsTagInfos") )

# prepare a path running the new modules
newJetTracksAssociator = cms.Sequence(
    newJetTracksAssociatorAtVertex
)

newJetBtaggingIP = cms.Sequence(
    newImpactParameterTagInfos * (
        newTrackCountingHighEffBJetTags +
        newTrackCountingHighPurBJetTags +
        newJetProbabilityBJetTags +
        newJetBProbabilityBJetTags
    )
)

newJetBtaggingSV = cms.Sequence(
    newImpactParameterTagInfos *
    newSecondaryVertexTagInfos * (
        newSimpleSecondaryVertexBJetTags +
        newCombinedSecondaryVertexBJetTags +
        newCombinedSecondaryVertexMVABJetTags
    )
)

newJetBtaggingEle = cms.Sequence(
    newSoftPFElectronsTagInfos *
    newSoftPFElectronBJetTags
)

newJetBtaggingMu = cms.Sequence(
    newSoftPFMuonsTagInfos *
    newSoftPFMuonBJetTags 
)

newJetBtagging = cms.Sequence(
    newJetBtaggingIP +
    newJetBtaggingSV +
    newJetBtaggingEle +
    newJetBtaggingMu
)
