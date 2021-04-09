// -*- C++ -*-
//
// Package:    L1Trigger/DemonstratorTools
// Class:      GTTFileReader
//
/**\class GTTFileReader GTTFileReader.cc L1Trigger/DemonstratorTools/plugins/GTTFileReader.cc

 Description: Example EDProducer class, illustrating how BoardDataReader can be used to
   read I/O buffer files (that have been created in hardware/firmware tests), decode
   the contained data, and store this in EDM collections.

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thomas Williams <thomas.williams@stfc.ac.uk>
//         Created:  Fri, 19 Feb 2021 01:10:55 GMT
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

#include "DataFormats/L1Trigger/interface/Vertex.h"
#include "L1Trigger/DemonstratorTools/interface/BoardDataReader.h"
#include "L1Trigger/DemonstratorTools/interface/codecs/vertices.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"

//
// class declaration
//

class GTTFileReader : public edm::stream::EDProducer<> {
public:
  explicit GTTFileReader(const edm::ParameterSet&);
  ~GTTFileReader() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  void endStream() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  l1t::demo::BoardDataReader fileReader_;
};

//
// constants, enums and typedefs
//

// NOTE: At least some of the info from these constants will eventually come from config files
constexpr size_t kGapLength(44);
constexpr size_t kVertexTMUX(6);
constexpr size_t kVertexChanIndex(0);

const std::map<size_t, l1t::demo::ChannelSpec> kChannelSpecs = {
    /* channel index -> {link TMUX, TMUX index, inter-packet gap} */
    {kVertexChanIndex, {kVertexTMUX, 0, kGapLength}}};

//
// static data member definitions
//

//
// constructors and destructor
//

GTTFileReader::GTTFileReader(const edm::ParameterSet& iConfig)
    : fileReader_(l1t::demo::parseFileFormat(iConfig.getUntrackedParameter<std::string>("format")),
                  iConfig.getParameter<std::vector<std::string>>("files"),
                  9,
                  6,
                  10,
                  kChannelSpecs) {
  produces<l1t::VertexWordCollection>();
}

GTTFileReader::~GTTFileReader() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void GTTFileReader::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  l1t::demo::BoardData eventData(fileReader_.getNextEvent());

  l1t::VertexWordCollection vertices(l1t::demo::codecs::decodeVertices(eventData.at(kVertexChanIndex)));

  std::cout << vertices.size() << " vertices found" << std::endl;

  iEvent.put(std::make_unique<l1t::VertexWordCollection>(vertices));

  /* this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  SetupData& setup = iSetup.getData(setupToken_);
*/
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void GTTFileReader::beginStream(edm::StreamID) {
  // please remove this method if not needed
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void GTTFileReader::endStream() {
  // please remove this method if not needed
}

// ------------ method called when starting to processes a run  ------------
/*
void
GTTFileReader::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
GTTFileReader::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
GTTFileReader::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
GTTFileReader::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GTTFileReader::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GTTFileReader);
