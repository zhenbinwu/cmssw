#include "FWCore/Framework/interface/MakerMacros.h"

#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <cmath>

#include "FWCore/Framework/interface/stream/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingMuon.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCalo.h"
#include "L1TriggerScouting/Utilities/interface/printScObjects.h"

using namespace l1ScoutingRun3;

// ----------------------------- CLASS DECLARATION  ----------------------------
class DumpScObjects : public edm::stream::EDAnalyzer<> {

  public:
    // constructor and destructor
    explicit DumpScObjects(const edm::ParameterSet&);
    ~DumpScObjects() override{};

    // method for analyzing the events
    void analyze(const edm::Event&, const edm::EventSetup&) override;

  private:

    // dump contenct of BX
    void printBx(unsigned bx);

    // the tokens to access the data
    edm::EDGetTokenT<ScMuonOrbitCollection>    gmtMuonsToken_;
    edm::EDGetTokenT<ScJetOrbitCollection>     caloJetsToken_;
    edm::EDGetTokenT<ScEGammaOrbitCollection>  caloEGammasToken_;
    edm::EDGetTokenT<ScTauOrbitCollection>     caloTausToken_;
    edm::EDGetTokenT<ScBxSumsOrbitCollection>   caloEtSumsToken_;

    edm::Handle<ScMuonOrbitCollection> muonHandle_;
    edm::Handle<ScJetOrbitCollection> jetHandle_;
    edm::Handle<ScEGammaOrbitCollection> eGammaHandle_;
    edm::Handle<ScTauOrbitCollection> tauHandle_;
    edm::Handle<ScBxSumsOrbitCollection> etSumHandle_;

    // the min and max BX to be analyzed
    unsigned minBx_;
    unsigned maxBx_;

    // select collection to be printed
    bool checkMuons_;
    bool checkJets_;
    bool checkEGammas_;
    bool checkTaus_;
    bool checkEtSums_;

    // dump a specific (ORBIT, BX RANGE)
    bool searchEvent_;
    unsigned orbitNum_;
    unsigned searchStartBx_;
    unsigned searchStopBx_;

    // utils
    bool skipEmptyBx_;
};
// -----------------------------------------------------------------------------


// -------------------------------- constructor  -------------------------------

DumpScObjects::DumpScObjects(const edm::ParameterSet& iConfig):
  minBx_(iConfig.getUntrackedParameter<unsigned>("minBx", 0)),
  maxBx_(iConfig.getUntrackedParameter<unsigned>("maxBx", 3564)),

  checkMuons_(iConfig.getUntrackedParameter<bool>("checkMuons", true)),
  checkJets_(iConfig.getUntrackedParameter<bool>("checkJets", true)),
  checkEGammas_(iConfig.getUntrackedParameter<bool>("checkEGammas", true)),
  checkTaus_(iConfig.getUntrackedParameter<bool>("checkTaus", true)),
  checkEtSums_(iConfig.getUntrackedParameter<bool>("checkEtSums", true)),

  searchEvent_(iConfig.getUntrackedParameter<bool>("searchEvent", false)),
  orbitNum_(iConfig.getUntrackedParameter<unsigned>("orbitNumber", 0)),
  searchStartBx_(iConfig.getUntrackedParameter<unsigned>("searchStartBx", 0)),
  searchStopBx_(iConfig.getUntrackedParameter<unsigned>("searchStopBx", 0)),

  skipEmptyBx_(iConfig.getUntrackedParameter<bool>("skipEmptyBx", true))
{

  if (checkMuons_) gmtMuonsToken_    = consumes<ScMuonOrbitCollection>(iConfig.getParameter<edm::InputTag>("gmtMuonsTag"));
  if (checkJets_) caloJetsToken_    = consumes<ScJetOrbitCollection>(iConfig.getParameter<edm::InputTag>("caloJetsTag"));
  if (checkEGammas_) caloEGammasToken_ = consumes<ScEGammaOrbitCollection>(iConfig.getParameter<edm::InputTag>("caloEGammasTag"));
  if (checkTaus_) caloTausToken_    = consumes<ScTauOrbitCollection>(iConfig.getParameter<edm::InputTag>("caloTausTag"));
  if (checkEtSums_) caloEtSumsToken_  = consumes<ScBxSumsOrbitCollection>(iConfig.getParameter<edm::InputTag>("caloEtSumsTag"));

}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void DumpScObjects::analyze(const edm::Event& iEvent, const edm::EventSetup& evSetup) {

  if (checkMuons_)   iEvent.getByToken(gmtMuonsToken_, muonHandle_);
  if (checkJets_)    iEvent.getByToken(caloJetsToken_, jetHandle_);
  if (checkEGammas_) iEvent.getByToken(caloEGammasToken_, eGammaHandle_);
  if (checkTaus_)    iEvent.getByToken(caloTausToken_, tauHandle_);
  if (checkEtSums_)  iEvent.getByToken(caloEtSumsToken_, etSumHandle_);

  // get the orbit number
  unsigned currOrbit = iEvent.id().event();

  // if we are looking for a specific orbit
  if (searchEvent_){
    if (currOrbit != orbitNum_) return;
    
    // found the orbit
    for (unsigned bx=searchStartBx_; bx<=searchStopBx_; bx++){
      printBx(bx);
    }
  } else {

    if (skipEmptyBx_){
      
      // create a set of non empty BXs
      std::set<unsigned> uniqueBx;

      if (checkMuons_) {
        for (const unsigned& bx: muonHandle_->getFilledBxs()){
          if ((bx>=minBx_) || (bx<=maxBx_)) uniqueBx.insert(bx);
        }
      }
      if (checkJets_) {
        for (const unsigned& bx: jetHandle_->getFilledBxs()){
          if ((bx>=minBx_) || (bx<=maxBx_)) uniqueBx.insert(bx);
        }
      }
      if (checkEGammas_) {
        for (const unsigned& bx: eGammaHandle_->getFilledBxs()){
          if ((bx>=minBx_) || (bx<=maxBx_)) uniqueBx.insert(bx);
        }
      }
      if (checkTaus_) {
        for (const unsigned& bx: tauHandle_->getFilledBxs()){
          if ((bx>=minBx_) || (bx<=maxBx_)) uniqueBx.insert(bx);
        }
      }
      if (checkEtSums_) {
        for (const unsigned& bx: etSumHandle_->getFilledBxs()){
          if ((bx>=minBx_) || (bx<=maxBx_)) uniqueBx.insert(bx);
        }
      }

      // process bx
      for (const unsigned& bx: uniqueBx){
        printBx(bx);
      }
      
    } 
    else {
      // dump all objects
      for (unsigned bx=minBx_; bx<=maxBx_; bx++){
        printBx(bx);
      }
    }
  }

}
// -----------------------------------------------------------------------------

void DumpScObjects::printBx(unsigned bx){
  std::cout << "BX = " << bx <<" ****" << std::endl;

  if(checkMuons_ && muonHandle_.isValid()){
    int i=0;
    for (auto muon = muonHandle_->begin(bx); muon!=muonHandle_->end(bx); muon++){
      std::cout  <<  "--- Muon " << i << " ---\n";
      printScMuon(*muon);
      i++;
    }
  }

  if(checkJets_ && jetHandle_.isValid()){
    int i=0;
    for (auto jet = jetHandle_->begin(bx); jet!=jetHandle_->end(bx); jet++){
      std::cout  <<  "--- Jet " << i << " ---\n";
      printScJet(*jet);
      i++;
    }
  }

  if(checkEGammas_ && jetHandle_.isValid()){
    int i=0;
    for (auto egamma = eGammaHandle_->begin(bx); egamma!=eGammaHandle_->end(bx); egamma++){
      std::cout  <<  "--- E/Gamma " << i << " ---\n";
      printScEGamma(*egamma);
      i++;
    }
  }

  if(checkTaus_ && tauHandle_.isValid()){
    int i=0;
    for (auto tau = tauHandle_->begin(bx); tau!=tauHandle_->end(bx); tau++){
      std::cout  <<  "--- Tau " << i << " ---\n";
      printScTau(*tau);
      i++;
    }
  }

  if(checkEtSums_ && etSumHandle_.isValid()){
    int i=0;
    for (auto sum = etSumHandle_->begin(bx); sum!=etSumHandle_->end(bx); sum++){
      std::cout  <<  "--- Calo Sums ---\n";
      printScBxSums(*sum);
      i++;
    }
  }

}

DEFINE_FWK_MODULE(DumpScObjects);