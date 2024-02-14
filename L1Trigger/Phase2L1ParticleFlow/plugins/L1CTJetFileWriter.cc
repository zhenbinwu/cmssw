#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/View.h"

#include "L1Trigger/DemonstratorTools/interface/BoardDataWriter.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"
#include "DataFormats/L1TParticleFlow/interface/PFJet.h"
#include "DataFormats/L1TParticleFlow/interface/gt_datatypes.h"

//
// class declaration
//

class L1CTJetFileWriter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit L1CTJetFileWriter(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ----------constants, enums and typedefs ---------
  std::vector<edm::ParameterSet> collections_;
  unsigned nJets_;
  size_t nFramesPerBX_;
  size_t ctl2BoardTMUX_;
  size_t gapLengthOutput_;
  size_t maxLinesPerFile_;
  std::map<l1t::demo::LinkId, std::pair<l1t::demo::ChannelSpec, std::vector<size_t>>> channelSpecsOutputToGT_;

  // ----------member functions ----------------------
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  std::vector<ap_uint<64>> encodeJets(const std::vector<l1t::PFJet> jets);

  l1t::demo::BoardDataWriter fileWriterOutputToGT_;
  std::vector<edm::EDGetTokenT<edm::View<l1t::PFJet>>> jetsTokens_;
};

L1CTJetFileWriter::L1CTJetFileWriter(const edm::ParameterSet& iConfig)
    : 
      collections_(iConfig.getUntrackedParameter<std::vector<edm::ParameterSet>>("collections",
                                                                                 std::vector<edm::ParameterSet>())),
      nJets_(iConfig.getParameter<unsigned>("nJets")),
      nFramesPerBX_(iConfig.getParameter<unsigned>("nFramesPerBX")),
      ctl2BoardTMUX_(iConfig.getParameter<unsigned>("TMUX")),
      gapLengthOutput_(ctl2BoardTMUX_ * nFramesPerBX_ - 2 * nJets_ * collections_.size()),
      maxLinesPerFile_(iConfig.getParameter<unsigned>("maxLinesPerFile")),
      channelSpecsOutputToGT_{{{"jets", 0}, {{ctl2BoardTMUX_, gapLengthOutput_}, {0}}}},
      fileWriterOutputToGT_(l1t::demo::parseFileFormat(iConfig.getParameter<std::string>("format")),
                            iConfig.getParameter<std::string>("outputFilename"),
                            iConfig.getParameter<std::string>("outputFileExtension"),
                            nFramesPerBX_,
                            ctl2BoardTMUX_,
                            maxLinesPerFile_,
                            channelSpecsOutputToGT_) {
        for(const auto &pset : collections_){
          jetsTokens_.push_back(consumes<edm::View<l1t::PFJet>>(pset.getParameter<edm::InputTag>("jets")));
        }
      }

void L1CTJetFileWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // 1) Pack collections in otder they're specified. jets then sums within collection
  std::vector<ap_uint<64>> link_words;
  for(const auto &token : jetsTokens_){
    // 2) Encode jet information onto vectors containing link data
    // TODO remove the sort here and sort the input collection where it's created
    const edm::View<l1t::PFJet>& jets = iEvent.get(token);
    std::vector<l1t::PFJet> sortedJets;
    sortedJets.reserve(jets.size());
    std::copy(jets.begin(), jets.end(), std::back_inserter(sortedJets));

    std::stable_sort(
        sortedJets.begin(), sortedJets.end(), [](l1t::PFJet i, l1t::PFJet j) { return (i.hwPt() > j.hwPt()); });
    const auto outputJets(encodeJets(sortedJets));
    link_words.insert(link_words.end(), outputJets.begin(), outputJets.end());
  }
  // 3) Pack jet information into 'event data' object, and pass that to file writer
  l1t::demo::EventData eventDataJets;
  eventDataJets.add({"jets", 0}, link_words);
  fileWriterOutputToGT_.addEvent(eventDataJets);
}

// ------------ method called once each job just after ending the event loop  ------------
void L1CTJetFileWriter::endJob() {
  // Writing pending events to file before exiting
  fileWriterOutputToGT_.flush();
}

std::vector<ap_uint<64>> L1CTJetFileWriter::encodeJets(const std::vector<l1t::PFJet> jets) {
  std::vector<ap_uint<64>> jet_words;
  for (unsigned i = 0; i < nJets_; i++) {
    l1t::PFJet j;
    if (i < jets.size()) {
      j = jets.at(i);
    } else {  // pad up to nJets_ with null jets
      l1t::PFJet j(0, 0, 0, 0, 0, 0);
    }
    jet_words.push_back(j.encodedJet()[0]);
    jet_words.push_back(j.encodedJet()[1]);
  }
  return jet_words;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1CTJetFileWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  {
    edm::ParameterSetDescription vpsd1;
    vpsd1.add<edm::InputTag>("jets", edm::InputTag("scPFL1PuppiEmulatorCorrected"));
    std::vector<edm::ParameterSet> temp1;
    temp1.reserve(1);
    {
      edm::ParameterSet temp2;
      temp2.addParameter<edm::InputTag>("jets", edm::InputTag("scPFL1PuppiEmulatorCorrected"));
      temp1.push_back(temp2);
    }
    desc.addVPSetUntracked("collections", vpsd1, temp1);
  }
  //desc.add<std::vector<edm::ParameterSet>>("collections");
  desc.add<std::string>("outputFilename");
  desc.add<std::string>("outputFileExtension", "txt");
  desc.add<uint32_t>("nJets", 12);
  desc.add<uint32_t>("nFramesPerBX", 9);
  desc.add<uint32_t>("TMUX", 6);
  desc.add<uint32_t>("maxLinesPerFile", 1024);
  desc.add<std::string>("format", "EMPv2");
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1CTJetFileWriter);
