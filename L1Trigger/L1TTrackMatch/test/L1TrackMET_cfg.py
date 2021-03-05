############################################################
# define basic process
############################################################

import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import os

############################################################
# edit options here
############################################################
L1TRK_INST ="L1TrackMET" ### if not in input DIGRAW then we make them in the above step
process = cms.Process(L1TRK_INST)

############################################################
# import standard configurations
############################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '')

############################################################
# input and output
############################################################

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))

readFiles = cms.untracked.vstring(
    "/store/relval/CMSSW_11_0_0/RelValTTbar_14TeV/GEN-SIM-DIGI-RAW/PU25ns_110X_mcRun4_realistic_v3_2026D49PU200-v2/10000/01054EE2-1B51-C449-91A2-5202A60D16A3.root"
)
secFiles = cms.untracked.vstring()

process.source = cms.Source ("PoolSource",
                            fileNames = readFiles,
                            secondaryFileNames = secFiles,
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            )


process.TFileService = cms.Service("TFileService", fileName = cms.string('CheckingJets_CMSSW11_CMS.root'), closeFileFast = cms.untracked.bool(True))

process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")
process.load("L1Trigger.L1TTrackMatch.L1TrackerEtMissProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1TrackerEtMissEmulatorProducer_cfi")
process.load("L1Trigger.L1TTrackMatch.L1TkMETAnalyser_cfi")

process.TTTracksEmulation = cms.Path(process.L1HybridTracks)
process.TTTracksEmulationWithTruth = cms.Path(process.L1HybridTracksWithAssociators)

process.pTkMET = cms.Path(process.L1TrackerEtMiss)
process.pTkEmuMET = cms.Path(process.L1TrackerEmuEtMiss)

############################################################
# Primary vertex
############################################################
process.load("L1Trigger.L1TTrackMatch.L1TkPrimaryVertexProducer_cfi")
process.pPV = cms.Path(process.L1TkPrimaryVertex)



process.analysis = cms.Path(process.L1TkMETAnalyser)

process.out = cms.OutputModule( "PoolOutputModule",
                                fastCloning = cms.untracked.bool( False ),
                                fileName = cms.untracked.string("test.root" )
		               )
process.pOut = cms.EndPath(process.out)

process.schedule = cms.Schedule(process.TTTracksEmulationWithTruth, process.pPV, process.pTkMET,process.pTkEmuMET, process.analysis)
