#ifndef PHASE2GMT_CONEVRTEDTTRACK
#define PHASE2GMT_CONEVRTEDTTRACK
#include "L1Trigger/Phase2L1GMT/interface/Constants.h"

namespace Phase2L1GMT {

  class ConvertedTTTrack {

  public:
  ConvertedTTTrack(const ap_uint<1>& charge,
		   const ap_int<BITSCURV>& curvature,
		   const ap_uint<BITSPT>& pt,
		   const ap_int<BITSETA>& eta,
		   const ap_int<BITSPHI>& phi,
		   const ap_int<BITSZ0>& z0,
		   const ap_int<BITSD0>& d0,
		   const ap_uint<BITSQUALITY>& quality): 
    charge_(charge),
    curvature_(curvature),  
    pt_(pt),
    eta_(eta),
    phi_(phi),  
    z0_(z0),
    d0_(d0),
    quality_(quality)  
    {
    }

    const ap_uint<1> charge() const {
      return charge_;
    }

    const ap_int<BITSCURV> curvature() const {
      return curvature_;
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
    const ap_int<BITSZ0> z0() const {
      return z0_;
    }
    const ap_int<BITSD0> d0() const {
      return d0_;
    }
    const ap_uint<BITSQUALITY> quality() const {
      return quality_;
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

    void print() const {
      printf("converted track charge=%d curvature=%d pt=%f,%d eta=%f,%d phi=%f,%d z0=%d d0=%d quality=%d\n",charge_.to_int(),curvature_.to_int(),offline_pt_,pt_.to_int(),offline_eta_,eta_.to_int(),offline_phi_,phi_.to_int(),z0_.to_int(),d0_.to_int(),quality_.to_int());
    }




  private:
    ap_uint<1>                charge_;
    ap_int<BITSCURV>          curvature_;
    ap_uint<BITSPT>           pt_;
    ap_int<BITSETA>           eta_;
    ap_int<BITSPHI>           phi_;
    ap_int<BITSZ0>            z0_;
    ap_int<BITSD0>            d0_;
    ap_uint<BITSQUALITY>      quality_;
    float                     offline_pt_;
    float                     offline_eta_;
    float                     offline_phi_;

  };
}

#endif
