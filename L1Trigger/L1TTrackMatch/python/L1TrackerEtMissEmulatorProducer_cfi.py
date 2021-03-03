import FWCore.ParameterSet.Config as cms

L1TrackerEmuEtMiss = cms.EDProducer('L1TrackerEtMissEmulatorProducer',
    L1TrackInputTag = cms.InputTag("TTTracksFromTrackletEmulation", "Level1TTTracks"),
    L1VertexInputTag = cms.InputTag("L1TkPrimaryVertex"),

    maxZ0 = cms.double ( 15. ) ,    # in cm
    maxEta = cms.double ( 2.4 ) ,   # max eta allowed for chosen tracks
    minPt = cms.double( 2. ),
    chi2rzdofMax = cms.double( 10. ), # max chi2/dof allowed for chosen tracks
    chi2rphidofMax = cms.double( 10. ), # max chi2/dof allowed for chosen tracks
    bendChi2Max = cms.double( 2. ),# max bendchi2 allowed for chosen tracks
    nStubsmin = cms.int32( 4 ),     # min number of stubs for the tracks

    nCordicSteps = cms.int32( 8 ), #Number of steps for cordic sqrt and phi computation
    Debug        = cms.int32( 5 ),  #0 - No Debug, 1 - LUT debug, 2 - Phi Debug, 3 - Z debug, 4 - Et Debug, 5 - Cordic Debug, 6 - Output
    WriteLUTs    = cms.bool(True)
    
)

