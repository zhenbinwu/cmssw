// -*- C++ -*-
//
// Package:    L1Trigger/L1TPrefireVetoBitFilter
// Class:      L1TPrefireVetoBitFilter
// 
/**\class L1TPrefireVetoBitFilter L1TPrefireVetoBitFilter.cc L1Trigger/L1TPrefireVetoBitFilter/plugins/L1TPrefireVetoBitFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Zhenbin Wu
//         Created:  Fri, 31 Jan 2020 21:14:55 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/L1TGlobal/interface/GlobalExtBlk.h"
//
// class declaration
//

class L1TPrefireVetoBitFilter : public edm::stream::EDFilter<> {
   public:
      explicit L1TPrefireVetoBitFilter(const edm::ParameterSet&);
      ~L1TPrefireVetoBitFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginStream(edm::StreamID) override;
      virtual bool filter(edm::Event&, const edm::EventSetup&) override;
      virtual void endStream() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  // ----------member data ---------------------------
  const edm::EDGetTokenT<GlobalAlgBlkBxCollection> token_;
  const edm::EDGetTokenT<GlobalExtBlkBxCollection> token_ext_;
  const unsigned int m_triggerRulePrefireVetoBit = 255;

      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1TPrefireVetoBitFilter::L1TPrefireVetoBitFilter(const edm::ParameterSet& iConfig)
      : token_(consumes<GlobalAlgBlkBxCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      token_ext_( consumes<GlobalExtBlkBxCollection>(iConfig.getParameter<edm::InputTag>("src_ext")))
{
   //now do what ever initialization is needed

}


L1TPrefireVetoBitFilter::~L1TPrefireVetoBitFilter()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
L1TPrefireVetoBitFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   //const std::vector<bool>* wordp = nullptr;
   bool unprefireable_bit = false;
   //edm::Handle<GlobalAlgBlkBxCollection> handleResults;
   //iEvent.getByToken(token_, handleResults);
   //wordp = &handleResults->at(0, 0).getAlgoDecisionFinal();
   edm::Handle<GlobalExtBlkBxCollection> handleExtResults;
   iEvent.getByToken(token_ext_, handleExtResults);
   unprefireable_bit = handleExtResults->at(0, 0).getExternalDecision(
           std::max(m_triggerRulePrefireVetoBit, GlobalExtBlk::maxExternalConditions - 1));


   return unprefireable_bit;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
L1TPrefireVetoBitFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void
L1TPrefireVetoBitFilter::endStream() {
}

// ------------ method called when starting to processes a run  ------------
/*
void
L1TPrefireVetoBitFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/
 
// ------------ method called when ending the processing of a run  ------------
/*
void
L1TPrefireVetoBitFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when starting to processes a luminosity block  ------------
/*
void
L1TPrefireVetoBitFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
L1TPrefireVetoBitFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1TPrefireVetoBitFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(L1TPrefireVetoBitFilter);
