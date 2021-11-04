import FWCore.ParameterSet.Config as cms

l1ctSeededConeJetFileWriter = cms.EDAnalyzer('L1CTJetFileWriter',
  jets = cms.InputTag("scPFL1PuppiEmulator"),
  outputFilename = cms.untracked.string("L1CTJetsPatterns"),
  format = cms.untracked.string("EMP")
)
