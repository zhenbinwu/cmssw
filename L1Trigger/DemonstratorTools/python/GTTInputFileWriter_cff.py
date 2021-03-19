import FWCore.ParameterSet.Config as cms

GTTInputFileWriter = cms.EDAnalyzer('GTTInputFileWriter',
  tracks = cms.untracked.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"),
  format = cms.untracked.string("APx")
)