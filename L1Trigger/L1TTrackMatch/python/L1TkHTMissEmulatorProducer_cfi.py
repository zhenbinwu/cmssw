import FWCore.ParameterSet.Config as cms

L1TrackerEmuHTMiss = cms.EDProducer("L1TkHTMissEmulatorProducer",
    L1TkJetInputTag = cms.InputTag("L1TrackJets", "L1TrackJets"),
    L1VertexInputTag = cms.InputTag("L1TkPrimaryVertex"),
    L1MHTCollectionName = cms.string("L1TrackerEmuHTMiss"),
    jet_maxEta = cms.double(2.4),
    jet_minPt = cms.double(5.0),
    jet_minNtracksLowPt=cms.int32(2),
    jet_minNtracksHighPt=cms.int32(3)
)
