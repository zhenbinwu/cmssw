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
  unsigned nJets;
  size_t kFramesPerBX; 
  size_t kCTL2BoardTMUX; 
  size_t kGapLengthOutput; 
  size_t kMaxLinesPerFile; 
  std::map<l1t::demo::LinkId, std::pair<l1t::demo::ChannelSpec, std::vector<size_t>>> kChannelSpecsOutputToGT; 

  // ----------member functions ----------------------
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  std::vector<ap_uint<64>> encodeJets(const edm::View<l1t::PFJet>& jets);

  edm::EDGetTokenT<edm::View<l1t::PFJet>> jetsToken_;
  l1t::demo::BoardDataWriter fileWriterOutputToGT_;
};

L1CTJetFileWriter::L1CTJetFileWriter(const edm::ParameterSet& iConfig) :
      nJets(iConfig.getParameter<unsigned>("nJets")),
      kFramesPerBX(iConfig.getParameter<unsigned>("nFramesPerBX")),
      kCTL2BoardTMUX(iConfig.getParameter<unsigned>("TMUX")),
      kGapLengthOutput(kCTL2BoardTMUX * kFramesPerBX - 2 * nJets),
      kMaxLinesPerFile(iConfig.getParameter<unsigned>("maxLinesPerFile")),
      kChannelSpecsOutputToGT{{{"jets", 0}, {{kCTL2BoardTMUX, kGapLengthOutput}, {0}}}},
      jetsToken_(consumes<edm::View<l1t::PFJet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      fileWriterOutputToGT_(l1t::demo::parseFileFormat(iConfig.getParameter<std::string>("format")),
                            iConfig.getParameter<std::string>("outputFilename"),
                            kFramesPerBX,
                            kCTL2BoardTMUX,
                            kMaxLinesPerFile,
                            kChannelSpecsOutputToGT) {
      }

void L1CTJetFileWriter::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // 1) Encode jet information onto vectors containing link data
  const auto outputJets(encodeJets(iEvent.get(jetsToken_)));

  // 2) Pack jet information into 'event data' object, and pass that to file writer
  l1t::demo::EventData eventDataJets;
  eventDataJets.add({"jets", 0}, outputJets);
  fileWriterOutputToGT_.addEvent(eventDataJets);

}

// ------------ method called once each job just after ending the event loop  ------------
void L1CTJetFileWriter::endJob() {
  // Writing pending events to file before exiting
  fileWriterOutputToGT_.flush();
}

std::vector<ap_uint<64>> L1CTJetFileWriter::encodeJets(const edm::View<l1t::PFJet>& jets) {
    std::vector<ap_uint<64>> jet_words;
    for(unsigned i=0; i < nJets; i++){
        l1t::PFJet j;
        if(i < jets.size()){
            j = jets.at(i);
        }else{ // pad up to nJets with null jets
            l1t::PFJet j(0,0,0,0,0,0);
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
  desc.add<edm::InputTag>("jets");
  desc.add<std::string>("outputFilename");
  desc.add<uint32_t>("nJets", 12);
  desc.add<uint32_t>("nFramesPerBX", 9);
  desc.add<uint32_t>("TMUX", 6);
  desc.add<uint32_t>("maxLinesPerFile", 1024);
  desc.add<std::string>("format","EMP");
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1CTJetFileWriter);
