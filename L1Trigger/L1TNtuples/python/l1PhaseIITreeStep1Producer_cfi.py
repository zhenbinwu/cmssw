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

   #GMT muons
   gmtMuonToken = cms.InputTag("L1SAMuonsGmt", "promptSAMuons"),  #("L1TkStubsGmt", ""), #("L1TkStubsGmt", ""),
   gmtTkMuonToken = cms.InputTag("L1TkMuonsGmt",""),


   scPFL1Puppi = cms.InputTag("scPFL1Puppi", ""),

   l1pfPhase1L1TJetToken  = cms.InputTag("Phase1L1TJetCalibrator9x9" ,   "Phase1L1TJetFromPfCandidates"), #use the 9x9 case
   l1pfPhase1L1TJetMET  = cms.InputTag("Phase1L1TJetProducer9x9" ,   "UncalibratedPhase1L1TJetFromPfCandidatesMET"), #use the 9x9 case
   l1pfPhase1L1TJetSums  = cms.InputTag("Phase1L1TJetSumsProducer9x9" ,   "Sums"), #use the 9x9 case

   caloJetToken = cms.InputTag("L1CaloJet","CaloJets"),
   caloJetHTTToken = cms.InputTag("L1CaloJetHTT","CaloJetHTT"),
   caloTauToken = cms.InputTag("L1CaloJet","CaloTaus"),
   L1HPSPFTauToken = cms.InputTag("HPSPFTauProducerPF",""),

   l1PFMet = cms.InputTag("L1MetPfProducer",""), #emulator

   zoPuppi = cms.InputTag("l1pfProducerBarrel","z0"),
   l1vertextdr = cms.InputTag("VertexProducer","l1vertextdr"),
   l1vertices = cms.InputTag("VertexProducer","l1vertices"),
   l1TkPrimaryVertex= cms.InputTag("L1VertexFinderEmulator","l1verticesEmulation"), #we need to rename this, but these are now emulated vertices!

   L1NNTauToken = cms.InputTag("l1NNTauProducerPuppi","L1PFTausNN"),
   L1NNTauPFToken = cms.InputTag("l1NNTauProducer","L1PFTausNN"),

   tkTrackerJetToken = cms.InputTag("L1TrackJetsEmulation", "L1TrackJets"),  #these are emulated
   tkTrackerJetDisplacedToken = cms.InputTag("L1TrackJetsExtendedEmulation", "L1TrackJetsExtended"), #emulated
	 
   tkMetToken = cms.InputTag("L1TrackerEmuEtMiss","L1TrackerEmuEtMiss"), #emulated
   tkMetDisplacedToken = cms.InputTag("L1TrackerEtMissExtended","L1TrackerExtendedEtMiss"),

   tkMhtTokens = cms.VInputTag( cms.InputTag("L1TrackerHTMiss","L1TrackerHTMiss")),
   tkMhtDisplacedTokens = cms.VInputTag( cms.InputTag("L1TrackerHTMissExtended","L1TrackerHTMissExtended")),

   maxL1Extra = cms.uint32(50)
)

#### Gen level tree

from L1Trigger.L1TNtuples.l1GeneratorTree_cfi  import l1GeneratorTree
genTree=l1GeneratorTree.clone()

runmenutree=cms.Path(l1PhaseIITree*genTree)




