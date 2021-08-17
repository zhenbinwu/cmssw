
/**\class L1TrackerHTMissEmulatorProducer L1TrackerHTMissEmulatorProducer.cc
 L1Trigger/L1TTrackMatch/plugins/L1TrackerHTMissEmulatorProducer.cc
 Description: Takes L1TTkJets and performs a integer emulation of Track-based missing HT, outputting a collection of EtSum 
*/

// system include files
#include <memory>
#include <numeric>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/L1TCorrelator/interface/TkPrimaryVertex.h"
#include "DataFormats/L1TCorrelator/interface/TkHTMiss.h"
#include "DataFormats/L1TCorrelator/interface/TkHTMissFwd.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

using namespace l1t;

class L1TkHTMissEmulatorProducer : public edm::stream::EDProducer<> {
public:
  explicit L1TkHTMissEmulatorProducer(const edm::ParameterSet&);
  ~L1TkHTMissEmulatorProducer() override;

private:
  virtual void beginJob();
  void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob(); 

  // ----------member data ---------------------------

  float jetMinPt_;            // [GeV]
  float jetMaxEta_;           // [rad]
  unsigned int minNtracksHighPt_;
  unsigned int minNtracksLowPt_;
  std::string L1MHTCollectionName_;

  const edm::EDGetTokenT<TkPrimaryVertexCollection> pvToken_;
  const edm::EDGetTokenT<TkJetCollection> jetToken_;
};

L1TkHTMissEmulatorProducer::L1TkHTMissEmulatorProducer(const edm::ParameterSet& iConfig)
  : pvToken_(consumes<TkPrimaryVertexCollection>(iConfig.getParameter<edm::InputTag>("L1VertexInputTag"))),
    jetToken_(consumes<TkJetCollection>(iConfig.getParameter<edm::InputTag>("L1TkJetInputTag"))) {
  jetMinPt_ = (float)iConfig.getParameter<double>("jet_minPt");
  jetMaxEta_ = (float)iConfig.getParameter<double>("jet_maxEta");
  minNtracksHighPt_ = iConfig.getParameter<int>("jet_minNtracksHighPt");
  minNtracksLowPt_ = iConfig.getParameter<int>("jet_minNtracksLowPt");
  // Name of output ED Product
  L1MHTCollectionName_ = (std::string)iConfig.getParameter<std::string>("L1MHTCollectionName");
  
  produces<std::vector<EtSum>>(L1MHTCollectionName_);

}

L1TkHTMissEmulatorProducer::~L1TkHTMissEmulatorProducer() {}

void L1TkHTMissEmulatorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  // std::unique_ptr<TkHTMissCollection> MHTCollection(new TkHTMissCollection);
  std::unique_ptr<std::vector<l1t::EtSum>> MHTCollection(new std::vector<l1t::EtSum>(0));

  // L1 primary vertex
  edm::Handle<TkPrimaryVertexCollection> L1VertexHandle;
  iEvent.getByToken(pvToken_, L1VertexHandle);
  std::vector<TkPrimaryVertex>::const_iterator vtxIter;

  // L1 track-trigger jets
  edm::Handle<TkJetCollection> L1TkJetsHandle;
  iEvent.getByToken(jetToken_, L1TkJetsHandle);
  std::vector<TkJet>::const_iterator jetIter;

  if (!L1TkJetsHandle.isValid()) {
    LogError("TkHTMissProducer") << "\nWarning: TkJetCollection not found in the event. Exit\n";
    return;
  }

  edm::Ref<TkPrimaryVertexCollection> L1VtxRef;  // null reference
  float sumPx = 0;
  float sumPy = 0;
  float HT = 0;

  // loop over jets
  for (jetIter = L1TkJetsHandle->begin(); jetIter != L1TkJetsHandle->end(); ++jetIter) {
    float tmp_jet_px = jetIter->px();
    float tmp_jet_py = jetIter->py();
    float tmp_jet_et = jetIter->et();
    float tmp_jet_pt = jetIter->pt();
    if (tmp_jet_pt < jetMinPt_)
      continue;
    if (fabs(jetIter->eta()) > jetMaxEta_)
      continue;
    if (jetIter->ntracks() < minNtracksLowPt_ && tmp_jet_et > 50)
      continue;
    if (jetIter->ntracks() < minNtracksHighPt_ && tmp_jet_et > 100)
      continue;
    sumPx += tmp_jet_px;
    sumPy += tmp_jet_py;
    HT += tmp_jet_pt;
  }  // end jet loop

  // define missing HT
  float et = sqrt(sumPx * sumPx + sumPy * sumPy);
  math::XYZTLorentzVector missingEt(-sumPx, -sumPy, 0, et);
  // edm::RefProd<TkJetCollection> jetCollRef(L1TkJetsHandle);
  // TkHTMiss tkHTM(missingEt, HT, jetCollRef, L1VtxRef);

  EtSum L1HTSum(missingEt, EtSum::EtSumType::kMissingHt, (int)et, 0, 0, 0);

  MHTCollection->push_back(L1HTSum);
  // iEvent.put(std::move(MHTCollection), "L1TrackerHTMiss");
  iEvent.put(std::move(MHTCollection), L1MHTCollectionName_);

} //end producer


void L1TkHTMissEmulatorProducer::beginJob() {}

void L1TkHTMissEmulatorProducer::endJob() {}

DEFINE_FWK_MODULE(L1TkHTMissEmulatorProducer);
