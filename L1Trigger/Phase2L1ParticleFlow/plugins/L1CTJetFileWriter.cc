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
#include "L1Trigger/Phase2L1ParticleFlow/src/newfirmware/dataformats/gt_datatypes.h"

//
// class declaration
//

class L1CTJetFileWriter : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit L1CTJetFileWriter(const edm::ParameterSet&);

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // ----------constants, enums and typedefs ---------
  // NOTE: At least some of the info from these constants will eventually come from config files
  static constexpr size_t kFramesPerTMUXPeriod = 9;
  static constexpr size_t kGapLengthOutput = 54-2*10;
  static constexpr size_t kCTL2BoardTMUX = 6;
  static constexpr size_t kMaxLinesPerFile = 1024;

  const std::map<l1t::demo::LinkId, std::pair<l1t::demo::ChannelSpec, std::vector<size_t>>>
      kChannelSpecsOutputToGT = {
          /* logical channel within time slice -> {{link TMUX, inter-packet gap}, vector of channel indices} */
          {{"jets", 0}, {{kCTL2BoardTMUX, kGapLengthOutput}, {0}}}};

  // ----------member functions ----------------------
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  std::array<std::vector<ap_uint<64>>, 1> encodeJets(const edm::View<l1t::PFJet>& jets);

  edm::EDGetTokenT<edm::View<l1t::PFJet>> jetsToken_;
  l1t::demo::BoardDataWriter fileWriterOutputToGT_;
};

L1CTJetFileWriter::L1CTJetFileWriter(const edm::ParameterSet& iConfig) :
      jetsToken_(consumes<edm::View<l1t::PFJet>>(iConfig.getUntrackedParameter<edm::InputTag>("jets"))),
      fileWriterOutputToGT_(l1t::demo::parseFileFormat(iConfig.getUntrackedParameter<std::string>("format")),
                            iConfig.getUntrackedParameter<std::string>("outputFilename"),
                            kFramesPerTMUXPeriod,
                            kCTL2BoardTMUX,
                            kMaxLinesPerFile,
                            kChannelSpecsOutputToGT) {}

void L1CTJetFileWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // 1) Encode jet information onto vectors containing link data
  const auto outputJets(encodeJets(iEvent.get(jetsToken_)));

  // 2) Pack jet information into 'event data' object, and pass that to file writer
  l1t::demo::EventData eventDataJets;
  eventDataJets.add({"jets", 0}, outputJets.at(0));
  fileWriterOutputToGT_.addEvent(eventDataJets);

}

// ------------ method called once each job just after ending the event loop  ------------
void L1CTJetFileWriter::endJob() {
  // Writing pending events to file before exiting
  fileWriterOutputToGT_.flush();
}

std::array<std::vector<ap_uint<64>>, 1> L1CTJetFileWriter::encodeJets(const edm::View<l1t::PFJet>& jets) {
    std::vector<ap_uint<64>> jet_words;
    std::for_each(jets.begin(), jets.end(), [&](auto jet){
        jet_words.push_back(jet.encodedJet()[0]);
        jet_words.push_back(jet.encodedJet()[1]);
    });
    std::array<std::vector<ap_uint<64>>, 1> link = {{ jet_words }};
    return link;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void L1CTJetFileWriter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1CTJetFileWriter);
