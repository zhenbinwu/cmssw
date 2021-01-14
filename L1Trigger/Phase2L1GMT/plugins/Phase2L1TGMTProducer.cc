#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "DataFormats/L1TCorrelator/interface/TkMuon.h"
#include "L1Trigger/Phase2L1GMT/interface/Processor.h"

//
// class declaration
//
using namespace Phase2L1GMT;


class Phase2L1TGMTProducer : public edm::stream::EDProducer<> {
   public:
      explicit Phase2L1TGMTProducer(const edm::ParameterSet&);
      ~Phase2L1TGMTProducer() override;

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      void beginStream(edm::StreamID) override;
      void produce(edm::Event&, const edm::EventSetup&) override;
      void endStream() override;
      std::unique_ptr<Processor> processor_;
      edm::EDGetTokenT<l1t::TkMuon::L1TTTrackCollection> srcTracks_;

  
  
  

};
Phase2L1TGMTProducer::Phase2L1TGMTProducer(const edm::ParameterSet& iConfig):
  processor_(new Processor(iConfig)),
  srcTracks_(consumes<l1t::TkMuon::L1TTTrackCollection>(iConfig.getParameter<edm::InputTag>("srcTracks")))
{
  produces <std::vector<l1t::TkMuon> >();
}


Phase2L1TGMTProducer::~Phase2L1TGMTProducer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}

 



//
// member functions
//

// ------------ method called to produce the data  ------------
void
Phase2L1TGMTProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   Handle<l1t::TkMuon::L1TTTrackCollection> trackHandle;
   iEvent.getByToken(srcTracks_,trackHandle);
   std::vector<edm::Ptr< l1t::TkMuon::L1TTTrackType > > tracks;
   for (uint i=0;i<trackHandle->size();++i) {
     edm::Ptr< l1t::TkMuon::L1TTTrackType > track(trackHandle, i);
     tracks.push_back(track);
   }

   std::vector<l1t::TkMuon> out = processor_->process(tracks);
   std::unique_ptr<std::vector<l1t::TkMuon> > out1 = std::make_unique<std::vector<l1t::TkMuon> >(out); 
   iEvent.put(std::move(out1));

}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
Phase2L1TGMTProducer::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
Phase2L1TGMTProducer::endStream() {
}

void
Phase2L1TGMTProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Phase2L1TGMTProducer);
