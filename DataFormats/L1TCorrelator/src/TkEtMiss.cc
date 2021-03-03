#include "DataFormats/L1TCorrelator/interface/TkEtMiss.h"

using namespace l1t;

TkEtMiss::TkEtMiss() {}

TkEtMiss::TkEtMiss(const LorentzVector& p4,
                   EtMissType type,
                   const double& etTotal,
                   const double& etMissPU,
                   const double& etTotalPU,
                   const edm::Ref<TkPrimaryVertexCollection>& avtxRef,
                   int bx)
    : L1Candidate(p4),
      type_(type),
      etTot_(etTotal),
      etMissPU_(etMissPU),
      etTotalPU_(etTotalPU),
      vtxRef_(avtxRef),
      bx_(bx) {}

TkEtMiss::TkEtMiss(const LorentzVector& p4,
                   EtMissType type,
                   const double& etTotal,
                   const double& etMissPU,
                   const double& etTotalPU,
                   int bx)
    : L1Candidate(p4), type_(type), etTot_(etTotal), etMissPU_(etMissPU), etTotalPU_(etTotalPU), bx_(bx) {}


TkEtMiss::TkEtMiss(const LorentzVector& p4,
                   EtMissType type,
                   const unsigned int& Etmiss,
                   const unsigned int& EtPhi,
                   const unsigned int& NumTracks,
                   const double& hwEtScale,
                   const double& hwPhiScale,
                   int bx)
    : L1Candidate(p4), type_(type), hwiEtmiss_(Etmiss), hwiEtphi_(EtPhi), hwNumTracks_(NumTracks) , hwEtScale_(hwEtScale),hwPhiScale_(hwPhiScale), bx_(bx) {}

TkEtMiss::TkEtMiss(const LorentzVector& p4,
                   EtMissType type,
                   const double&  Etmiss,
                   const double&  EtPhi,
                   const unsigned int& NumTracks ,
                   int bx)
    : L1Candidate(p4), type_(type), hwEtmiss_(Etmiss), hwEtphi_(EtPhi), hwNumTracks_(NumTracks) , bx_(bx) {}


void TkEtMiss::ihwtofloat(){
  hwEtmiss_ = hwEtScale_*hwiEtmiss_;
  hwEtphi_  = hwPhiScale_*hwiEtphi_;
}

void TkEtMiss::hwfloattoi(){
  hwiEtmiss_ = hwEtmiss_/ hwEtScale_;
  hwiEtphi_  = hwEtphi_ / hwPhiScale_;
}