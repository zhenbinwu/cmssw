import FWCore.ParameterSet.Config as cms

l1PhaseIITree = cms.EDAnalyzer("L1PhaseIITreeStep1Producer",

   egTokenBarrel = cms.InputTag("L1EGammaClusterEmuProducer",""),  #keep as is, not produced by GCT
   tkEGTokenBarrel = cms.InputTag("l1ctLayer1EG","L1TkEleEB"),
   tkEMTokenBarrel = cms.InputTag("l1ctLayer1EG","L1TkEmEB"),

   egTokenHGC = cms.InputTag("l1ctLayer1EG","L1EgEE"),
   tkEGTokenHGC = cms.InputTag("l1ctLayer1EG","L1TkEleEE"),
   tkEMTokenHGC = cms.InputTag("l1ctLayer1EG","L1TkEmEE"),


   muonKalman = cms.InputTag("simKBmtfDigis","BMTF"),
   muonOverlap = cms.InputTag("simOmtfDigis","OMTF"),
   muonEndcap = cms.InputTag("simEmtfDigis",""),
   TkMuonToken = cms.InputTag("L1TkMuons",""),

   #Global muons
   muonToken = cms.untracked.InputTag("simGmtStage2Digis", ""),
   TkGlbMuonToken = cms.InputTag("L1TkGlbMuons",""),

   scPFL1Puppi = cms.InputTag("scPFL1Puppi", ""),

   l1pfPhase1L1TJetToken  = cms.InputTag("Phase1L1TJetCalibrator" ,   "Phase1L1TJetFromPfCandidates"), # not there yet either

   caloTauToken = cms.InputTag("L1CaloJet","CaloTaus"),
   L1HPSPFTauToken = cms.InputTag("HPSPFTauProducerPF",""),

   l1PFMet = cms.InputTag("L1MetPfProducer",""),

   zoPuppi = cms.InputTag("l1pfProducerBarrel","z0"),
   l1vertextdr = cms.InputTag("VertexProducer","l1vertextdr"),
   l1vertices = cms.InputTag("VertexProducer","l1vertices"),
   l1TkPrimaryVertex= cms.InputTag("L1VertexFinderEmulator","l1verticesEmulation"), #we need to rename this, but these are now emulated vertices!

   L1NNTauToken = cms.InputTag("l1NNTauProducerPuppi","L1PFTausNN"),
   L1NNTauPFToken = cms.InputTag("l1NNTauProducer","L1PFTausNN"),

   tkTrackerJetToken = cms.InputTag("L1TrackJetsEmulated", "L1TrackJets"),
   tkTrackerJetDisplacedToken = cms.InputTag("L1TrackJetsEmulatedExtended", "L1TrackJetsExtended"),

   tkMetToken = cms.InputTag("L1TrackerEtMiss","L1TrackerEtMiss"),
   tkMetDisplacedToken = cms.InputTag("L1TrackerEtMissExtended","L1TrackerEtMissExtended"),
   
   tkMhtTokens = cms.VInputTag( cms.InputTag("L1TrackerHTMiss","L1TrackerHTMiss")),
   tkMhtDisplacedTokens = cms.VInputTag( cms.InputTag("L1TrackerHTMissExtended","L1TrackerHTMissExtended")),

   maxL1Extra = cms.uint32(50)
)

#### Gen level tree

from L1Trigger.L1TNtuples.l1GeneratorTree_cfi  import l1GeneratorTree
genTree=l1GeneratorTree.clone()

runmenutree=cms.Path(l1PhaseIITree*genTree)




