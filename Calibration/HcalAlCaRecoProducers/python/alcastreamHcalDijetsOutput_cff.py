import FWCore.ParameterSet.Config as cms

# output block for alcastream HCAL Dijets
# output module 
#  module alcastreamHcalDijetsOutput = PoolOutputModule
alcastreamHcalDijetsOutput = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_DiJProd_*_*')
)

