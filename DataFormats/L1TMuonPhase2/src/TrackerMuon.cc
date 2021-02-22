
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"

using namespace l1t;

TrackerMuon::TrackerMuon():
  word_(0) { }


TrackerMuon::TrackerMuon(const edm::Ptr<L1TTTrackType>& trk,
			 bool charge,
			 uint pt,
			 int eta,
			 int phi,
			 int z0,
			 int d0,
			 uint quality):
  L1Candidate(LorentzVector(trk->momentum().x(),trk->momentum().y(),trk->momentum().z(),trk->momentum().mag()),pt,eta,phi,quality),
  trkPtr_(trk),
  hwCharge_(charge),
  hwZ0_(z0),
  hwD0_(d0),
  hwBeta_(15),
  word_(0) {


}


TrackerMuon::~TrackerMuon() {}


void TrackerMuon::print() const {
  printf("Tracker Muon: charge=%d pt=%d,%f eta=%d,%f phi=%d,%f z0=%d d0=%d isolation=%d beta=%d quality=%d\n",hwCharge_,hwPt(),p4().pt(),hwEta(),p4().eta(),hwPhi(),p4().phi(),hwZ0_,hwD0_,hwIso(),hwBeta_,hwQual());
}






  


