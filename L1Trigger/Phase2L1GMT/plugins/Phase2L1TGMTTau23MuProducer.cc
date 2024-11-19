// -*- C++ -*-
//
// Package:    L1Trigger/Phase2L1GMT
// Class:      Phase2L1TGMTTau23MuProducer
//
/**\class Phase2L1TGMTTau23MuProducer Phase2L1TGMTTau23MuProducer.cc L1Trigger/Phase2L1GMT/plugins/Phase2L1TGMTTau23MuProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhenbin Wu
//         Created:  Wed, 05 Jun 2024 21:46:26 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/L1TMuonPhase2/interface/Tau23Mu.h"
#include "L1Trigger/Phase2L1GMT/interface/Tauto3Mu.h"

//
// class declaration
//
using namespace Phase2L1GMT;

class Phase2L1TGMTTau23MuProducer : public edm::stream::EDProducer<> {
public:
  explicit Phase2L1TGMTTau23MuProducer(const edm::ParameterSet&);
  ~Phase2L1TGMTTau23MuProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;
  std::unique_ptr<Tauto3Mu> tau_;
  edm::EDGetTokenT<l1t::TrackerMuonCollection> srcTracks_;
  int minTrackStubs_;
};

Phase2L1TGMTTau23MuProducer::Phase2L1TGMTTau23MuProducer(const edm::ParameterSet& iConfig) 
  : tau_(new Tauto3Mu(iConfig)),
  srcTracks_(consumes<std::vector<l1t::TrackerMuon> >(iConfig.getParameter<edm::InputTag>("srcMuons")))
  //srcStubs_(consumes<std::vector<l1t::MuonStub> >(iConfig.getParameter<edm::InputTag>("srcStubs"))),
{
  produces<std::vector<l1t::Tau23Mu> >();
}

Phase2L1TGMTTau23MuProducer::~Phase2L1TGMTTau23MuProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void Phase2L1TGMTTau23MuProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  Handle<l1t::TrackerMuonCollection> trackHandle;
  iEvent.getByToken(srcTracks_, trackHandle);

  std::vector<l1t::Tau23Mu> out = tau_->GetTau3Mu(*trackHandle);
  std::unique_ptr<std::vector<l1t::Tau23Mu> > out1 = std::make_unique<std::vector<l1t::Tau23Mu> >(out);
  iEvent.put(std::move(out1));
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void Phase2L1TGMTTau23MuProducer::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void Phase2L1TGMTTau23MuProducer::endStream() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Phase2L1TGMTTau23MuProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1TGMTTau23MuProducer);
