// -*- C++ -*-
//
// Package:    L1Trigger/DemonstratorTools
// Class:      GTTInputFileWriter
//
/**\class GTTInputFileWriter GTTInputFileWriter.cc L1Trigger/DemonstratorTools/plugins/GTTInputFileWriter.cc

 Description: Example EDAnalyzer class, illustrating how BoardDataWriter can be used to write I/O buffer files for hardware/firmware tests

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Thomas Williams <thomas.williams@stfc.ac.uk>
//         Created:  Mon, 15 Feb 2021 00:39:44 GMT
//
//

// system include files
#include <memory>

#include "ap_int.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/Common/interface/View.h"

#include "L1Trigger/DemonstratorTools/interface/BoardDataWriter.h"
#include "L1Trigger/DemonstratorTools/interface/codecs/tracks.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

class GTTInputFileWriter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit GTTInputFileWriter(const edm::ParameterSet&);
  ~GTTInputFileWriter() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  typedef TTTrack<Ref_Phase2TrackerDigi_> Track_t;

  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  edm::EDGetTokenT<edm::View<Track_t>> tracksToken_;
  size_t eventCount_;
  l1t::demo::BoardDataWriter fileWriter_;
};

//
// constants, enums and typedefs
//

// NOTE: At least some of the info from these constants will eventually come from config files
constexpr size_t kGapLength(6);
constexpr size_t kTrackTMUX(18);

const std::map<size_t, l1t::demo::ChannelSpec> kChannelSpecs = {
    /* channel index -> {link TMUX, TMUX index, inter-packet gap} */
    {0, {kTrackTMUX, 0, kGapLength}},   {1, {kTrackTMUX, 0, kGapLength}},   {2, {kTrackTMUX, 0, kGapLength}},
    {3, {kTrackTMUX, 0, kGapLength}},   {4, {kTrackTMUX, 0, kGapLength}},   {5, {kTrackTMUX, 0, kGapLength}},
    {6, {kTrackTMUX, 0, kGapLength}},   {7, {kTrackTMUX, 0, kGapLength}},   {8, {kTrackTMUX, 0, kGapLength}},

    {9, {kTrackTMUX, 0, kGapLength}},   {10, {kTrackTMUX, 0, kGapLength}},  {11, {kTrackTMUX, 0, kGapLength}},
    {12, {kTrackTMUX, 0, kGapLength}},  {13, {kTrackTMUX, 0, kGapLength}},  {14, {kTrackTMUX, 0, kGapLength}},
    {15, {kTrackTMUX, 0, kGapLength}},  {16, {kTrackTMUX, 0, kGapLength}},  {17, {kTrackTMUX, 0, kGapLength}},

    {18, {kTrackTMUX, 6, kGapLength}},  {19, {kTrackTMUX, 6, kGapLength}},  {20, {kTrackTMUX, 6, kGapLength}},
    {21, {kTrackTMUX, 6, kGapLength}},  {22, {kTrackTMUX, 6, kGapLength}},  {23, {kTrackTMUX, 6, kGapLength}},
    {24, {kTrackTMUX, 6, kGapLength}},  {25, {kTrackTMUX, 6, kGapLength}},  {26, {kTrackTMUX, 6, kGapLength}},

    {27, {kTrackTMUX, 6, kGapLength}},  {28, {kTrackTMUX, 6, kGapLength}},  {29, {kTrackTMUX, 6, kGapLength}},
    {30, {kTrackTMUX, 6, kGapLength}},  {31, {kTrackTMUX, 6, kGapLength}},  {32, {kTrackTMUX, 6, kGapLength}},
    {33, {kTrackTMUX, 6, kGapLength}},  {34, {kTrackTMUX, 6, kGapLength}},  {35, {kTrackTMUX, 6, kGapLength}},

    {36, {kTrackTMUX, 12, kGapLength}}, {37, {kTrackTMUX, 12, kGapLength}}, {38, {kTrackTMUX, 12, kGapLength}},
    {39, {kTrackTMUX, 12, kGapLength}}, {40, {kTrackTMUX, 12, kGapLength}}, {41, {kTrackTMUX, 12, kGapLength}},
    {42, {kTrackTMUX, 12, kGapLength}}, {43, {kTrackTMUX, 12, kGapLength}}, {44, {kTrackTMUX, 12, kGapLength}},

    {45, {kTrackTMUX, 12, kGapLength}}, {46, {kTrackTMUX, 12, kGapLength}}, {47, {kTrackTMUX, 12, kGapLength}},
    {48, {kTrackTMUX, 12, kGapLength}}, {49, {kTrackTMUX, 12, kGapLength}}, {50, {kTrackTMUX, 12, kGapLength}},
    {51, {kTrackTMUX, 12, kGapLength}}, {52, {kTrackTMUX, 12, kGapLength}}, {53, {kTrackTMUX, 12, kGapLength}}};

//
// static data member definitions
//

//
// constructors and destructor
//

GTTInputFileWriter::GTTInputFileWriter(const edm::ParameterSet& iConfig)
    : tracksToken_(consumes<edm::View<Track_t>>(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))),
      eventCount_(0),
      fileWriter_(l1t::demo::parseFileFormat(iConfig.getUntrackedParameter<std::string>("format")),
                  iConfig.getUntrackedParameter<std::string>("outputFilename"),
                  9,
                  6,
                  1024,
                  kChannelSpecs) {
  //now do what ever initialization is needed
}

GTTInputFileWriter::~GTTInputFileWriter() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

void GTTInputFileWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Select links for correct time slice (18 input links per time slice)
  // TODO: Find nicer solution - e.g. having BoardDataWriter automatically shuffle
  //       data between channels based on event-independent 'logical' ID
  //       (at least nicer solution than storing event count as member var)
  const size_t baseIdx = (eventCount_ % 3) * 18;
  eventCount_++;

  // 1) Encode track information onto vectors containing link data
  const std::array<l1t::demo::BoardData::Channel, 18> trackData(
      l1t::demo::codecs::encodeTracks(iEvent.get(tracksToken_)));

  // 2) Pack track information into 'board data' object, and pass that to file writer
  l1t::demo::BoardData boardData;
  for (size_t i = 0; i < 18; i++)
    boardData.add(baseIdx + i, trackData.at(i));

  fileWriter_.addEvent(boardData);
}

// ------------ method called once each job just before starting event loop  ------------
void GTTInputFileWriter::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void GTTInputFileWriter::endJob() {
  // Writing pending events to file before exiting
  fileWriter_.flush();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GTTInputFileWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GTTInputFileWriter);
