#ifndef PHASE2GMT_MUONROI
#define PHASE2GMT_MUONROI
#include "L1Trigger/Phase2L1GMT/interface/Constants.h"
#include "DataFormats/L1TMuonPhase2/interface/L1TPhase2GMTStub.h"

namespace Phase2L1GMT {

  class MuonROI {

  public:
  MuonROI(const ap_uint<1>& charge,
	  const ap_uint<BITSPT>& pt,
	  const ap_int<BITSETA>& eta,
	  const ap_int<BITSPHI>& phi,
	  const ap_uint<BITSMUONQUALITY>& quality,
	  const ap_uint<3>& station
	  ): 
    charge_(charge),
    pt_(pt),
    eta_(eta),
    phi_(phi),  
    quality_(quality),  
    referenceStation_(station)
    {
    }

    const ap_uint<1> charge() const {
      return charge_;
    }

    const ap_uint<BITSPT> pt() const {
      return pt_;
    }

    const ap_int<BITSETA> eta() const {
      return eta_;
    }
    const ap_int<BITSPHI> phi() const {
      return phi_;
    }
    const ap_int<BITSMUONQUALITY> quality() const {
      return quality_;
    }
    const ap_uint<3> referenceStation() const {
      return referenceStation_;
    }
    }
    const float offline_pt() const {
      return offline_pt_;
    }
    const float offline_eta() const {
      return offline_eta_;
    }
    const float offline_phi() const {
      return offline_phi_;
    }


    void setOfflineQuantities(float pt,float eta, float phi) {
      offline_pt_=pt;
      offline_eta_=eta;
      offline_phi_=phi;
    }

    void addStub(const L1TPhase2GMTStubRef& stub){
      stubs_.push_back(stub);
    }

    const L1TPhase2GMTStubRefVector& stubs() {
      return stubs_;
    } 
  private:
    ap_uint<1>                charge_;
    ap_uint<BITSPT>           pt_;
    ap_int<BITSETA>           eta_;
    ap_int<BITSPHI>           phi_;
    ap_uint<BITSMUONQUALITY>      quality_;
    ap_uint<3>                referenceStation_;
    float                     offline_pt_;
    float                     offline_eta_;
    float                     offline_phi_;
    L1TPhase2GMTStubRefVector stubs_;

  };
}

#endif
