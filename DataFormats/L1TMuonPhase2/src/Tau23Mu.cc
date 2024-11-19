// ===========================================================================
// 
//       Filename:  Tau23Mu.cc
// 
//    Description: 
// 
//        Version:  1.0
//        Created:  09/22/2021 04:53:29 PM
//       Compiler:  g++ -std=c++11
// 
//         Author:  Zhenbin Wu (benwu)
//          Email:  zhenbin.wu@gmail.com
// 
// ===========================================================================

#include "DataFormats/L1TMuonPhase2/interface/Tau23Mu.h"

using namespace l1t;

Tau23Mu::Tau23Mu() : word_(0) {}

Tau23Mu::Tau23Mu(uint pt, int eta, int phi, uint mass, uint quality, int mu1, int mu2, int mu3)
    : L1Candidate(math::PtEtaPhiMLorentzVector {0, 0, 0, 0}, pt, eta, phi, quality),
    hwMass_(mass), mu1_(mu1), mu2_(mu2), mu3_(mu3), word_(0){}

Tau23Mu::~Tau23Mu() {}

void Tau23Mu::print() const {
  LogDebug("Tau23Mu") << "Tau to 3 muon: pt=" << hwPt() << "," << p4().pt()
                     << " eta=" << hwEta() << "," << p4().eta() << " phi=" << hwPhi() << "," << p4().phi()
                     << " mass=" << hwMass_;
}
