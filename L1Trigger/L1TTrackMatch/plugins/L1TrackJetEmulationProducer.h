#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <cstdlib>
#include "DataFormats/L1Trigger/interface/TkJetWord.h"

using namespace std;

namespace ExtraBits {
  const int PT_BITS = 18;
  const int ETA_BITS = 20;
  const int PHI_BITS = 21;
  const int Z0_BITS = 22;
}

typedef ap_ufixed<14 + ExtraBits::PT_BITS, 12, AP_TRN, AP_SAT> pt_intern;
typedef ap_int<12+ExtraBits::ETA_BITS> glbeta_intern;
typedef ap_int<11+ExtraBits::PHI_BITS> glbphi_intern;
typedef ap_int<10+ExtraBits::Z0_BITS> z0_intern;         // 40cm / 0.1

namespace Convert {
  const int INTPHI_PI = 720;
  const int INTPHI_TWOPI = 2 * INTPHI_PI;
  constexpr float INTPT_LSB_POW = pow(2.0,-2 - ExtraBits::PT_BITS);
  constexpr float INTPT_LSB = INTPT_LSB_POW;
  constexpr float ETA_LSB_POW = pow(2.0,-1 * ExtraBits::ETA_BITS);
  constexpr float ETA_LSB = M_PI / INTPHI_PI * ETA_LSB_POW;
  constexpr float PHI_LSB_POW = pow(2.0,-1 * ExtraBits::PHI_BITS);
  constexpr float PHI_LSB = M_PI / INTPHI_PI * PHI_LSB_POW;
  constexpr float Z0_LSB_POW = pow(2.0,-1 * ExtraBits::Z0_BITS);
  constexpr float Z0_LSB = 0.05 * Z0_LSB_POW;
  inline float floatPt(pt_intern pt) { return pt.to_float(); }
  inline int intPt(pt_intern pt) { return (ap_ufixed<16+ExtraBits::PT_BITS, 14>(pt)).to_int(); }
  inline float floatEta(glbeta_intern eta) { return eta.to_float() * ETA_LSB; }
  inline float floatPhi(glbphi_intern phi) { return phi.to_float() * PHI_LSB; }
  inline float floatZ0(z0_intern z0) { return z0.to_float() * Z0_LSB; }

  inline pt_intern makePt(int pt) { return ap_ufixed<16+ExtraBits::PT_BITS, 14>(pt); }
  inline pt_intern makePtFromFloat(float pt) { return pt_intern(INTPT_LSB_POW * round(pt / INTPT_LSB_POW)); }
  inline z0_intern makeZ0(float z0) { return z0_intern(round(z0 / Z0_LSB)); }

  inline ap_uint<pt_intern::width> ptToInt(pt_intern pt) {
    // note: this can be synthethized, e.g. when pT is used as intex in a LUT
    ap_uint<pt_intern::width> ret = 0;
    ret(pt_intern::width - 1, 0) = pt(pt_intern::width - 1, 0);
    return ret;
  }

  inline glbeta_intern makeGlbEta(float eta) { return round(eta / ETA_LSB); }
  inline glbeta_intern makeGlbEtaRoundEven(float eta) {
    glbeta_intern ghweta = round(eta / ETA_LSB);
    return (ghweta % 2) ? glbeta_intern(ghweta + 1) : ghweta;
  }

  inline glbphi_intern makeGlbPhi(float phi) { return round(phi / PHI_LSB); }

};  // namespace Scales


//Each individual box in the eta and phi dimension.
//  Also used to store final cluster data for each zbin.
struct EtaPhiBin {
  pt_intern pTtot;
  l1t::nt_t ntracks;
  l1t::nx_t nxtracks;
  bool used;
  glbphi_intern phi;  //average phi value (halfway b/t min and max)
  glbeta_intern eta;  //average eta value
};

//store important information for plots
struct MaxZBin {
  int znum;    //Numbered from 0 to nzbins (16, 32, or 64) in order
  int nclust;  //number of clusters in this bin
  z0_intern zbincenter;
  EtaPhiBin *clusters;  //list of all the clusters in this bin
  pt_intern ht;         //sum of all cluster pTs--only the zbin with the maximum ht is stored
};
