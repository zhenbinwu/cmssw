import FWCore.ParameterSet.Config as cms

GTTFileWriter = cms.EDAnalyzer('GTTFileWriter',
  tracks = cms.untracked.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"),
  convertedTracks = cms.untracked.InputTag("L1GTTInputProducer", "Level1TTTracksConverted"),
  vertices = cms.untracked.InputTag("VertexProducer", "l1verticesEmulation"),
  jets = cms.untracked.InputTag("L1TrackJetsEmulation","L1TrackJets"),
  htmiss = cms.untracked.InputTag("L1TrackerEmuHTMiss", "L1TrackerEmuHTMiss"),
  etmiss = cms.untracked.InputTag("L1TrackerEmuEtMiss", "L1TrackerEmuEtMiss"),
  inputFilename = cms.untracked.string("L1GTTInputFile"),
  inputConvertedFilename = cms.untracked.string("L1GTTInputConvertedFile"),
  outputCorrelatorFilename = cms.untracked.string("L1GTTOutputToCorrelatorFile"),
  outputGlobalTriggerFilename = cms.untracked.string("L1GTTOutputToGlobalTriggerFile"),
  format = cms.untracked.string("APx")
)