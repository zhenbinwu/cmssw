import FWCore.ParameterSet.Config as cms

GTTOutputFileReader = cms.EDProducer('GTTOutputFileReader',
  files = cms.vstring("gttOutput_0.txt") #, "gttOutput_1.txt")
)