
import FWCore.ParameterSet.Config as cms
import FWCore.Utilities.FileUtils as FileUtils
import FWCore.ParameterSet.VarParsing as VarParsing


# PART 1 : PARSE ARGUMENTS

options = VarParsing.VarParsing ('analysis')
options.parseArguments()

inputFiles = []
for filePath in options.inputFiles:
    if filePath.endswith(".root"):
        inputFiles.append(filePath)
    else:
        inputFiles += FileUtils.loadListFromFile(filePath)


# PART 2: SETUP MAIN CMSSW PROCESS 

process = cms.Process("GTTOutputValidation")

process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D49_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(inputFiles) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

process.load("L1Trigger.TrackFindingTracklet.L1HybridEmulationTracks_cff")
process.load('L1Trigger.DemonstratorTools.GTTOutputFileReader_cff')
process.GTTOutputFileReader.files = cms.vstring("test/gtt/example_vertex_apx.txt")
process.GTTOutputFileReader.format = cms.untracked.string("apx")

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.Timing = cms.Service("Timing", summaryOnly = cms.untracked.bool(True))

process.p = cms.Path(process.L1HybridTracks * process.GTTOutputFileReader) # vertex emulator & FW-emulator comparsion module need to be added here
