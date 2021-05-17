#include "DataFormats/L1TMuonPhase2/interface/SAMuon.h"

using namespace l1t;

SAMuon::SAMuon() : hwZ0_(0), hwD0_(0), word_(0) {}

SAMuon::SAMuon(const l1t::Muon& mu, bool charge, uint pt, int eta, int phi, int z0, int d0, uint quality)
    : L1Candidate(mu.p4(), pt, eta, phi, quality), hwCharge_(charge), hwZ0_(z0), hwD0_(d0), word_(0) {}

SAMuon::~SAMuon() {}

void SAMuon::print() const {
  printf("Standalone Muon: charge=%d pt=%d,%f eta=%d,%f phi=%d,%f z0=%d d0=%d isolation=%d beta=%d quality=%d\n",
         hwCharge_,
         hwPt(),
         p4().pt(),
         hwEta(),
         p4().eta(),
         hwPhi(),
         p4().phi(),
         hwZ0_,
         hwD0_,
         hwIso(),
         hwBeta_,
         hwQual());
}
