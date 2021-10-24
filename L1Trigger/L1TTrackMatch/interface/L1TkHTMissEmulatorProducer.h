#ifndef L1Trigger_L1TTrackMatch_L1TkHTMissEmulatorProducer_HH
#define L1Trigger_L1TTrackMatch_L1TkHTMissEmulatorProducer_HH 

// Original Author:  Hardik Routray                                                                                               
//         Created:  Mon, 11 Oct 2021 

#include <ap_int.h>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1Trigger/interface/TkJetWord.h"

using namespace std;
// Namespace that defines constants and types used by the HTMiss Emulation

namespace l1tmhtemu { 

  const unsigned int kPtSize_ = l1t::TkJetWord::TkJetBitWidths::kPtSize;
  const unsigned int kEtaSize_ = l1t::TkJetWord::TkJetBitWidths::kGlbEtaSize;
  const unsigned int kPhiSize_ = l1t::TkJetWord::TkJetBitWidths::kGlbPhiSize;

  const unsigned int kInternalPtWidth{kPtSize_};
  const unsigned int kInternalEtaWidth{kEtaSize_};
  const unsigned int kInternalPhiWidth{kPhiSize_};

  // extra room for sumPx, sumPy
  const unsigned int kEtExtra{10};

  const unsigned int kMHTSize{15};  // For output Magnitude default 15
  const unsigned int kMHTPhiSize{14}; // For output Phi default 14
  const float kMaxMHT{4096};  // 4 TeV
  const float kMaxMHTPhi{2 * M_PI};

  typedef ap_uint<5> ntracks_t;
  typedef ap_uint<kInternalPtWidth> pt_t;
  typedef ap_int<kInternalEtaWidth> eta_t;
  typedef ap_int<kInternalPhiWidth> phi_t;

  typedef ap_int<kInternalPtWidth + kEtExtra> Et_t;
  typedef ap_uint<kMHTSize> MHT_t;
  typedef ap_uint<kMHTPhiSize> MHTphi_t;

  const unsigned int kMHTBins = 1 << kMHTSize;
  const unsigned int kMHTPhiBins = 1 << kMHTPhiSize;

  const double kStepPt_ = 0.25;
  const double kStepEta_ = M_PI / (720);
  const double kStepPhi_ = M_PI / (720);

  const double kStepPt{kStepPt_}; 
  const double kStepEta{kStepEta_};
  const double kStepPhi{kStepPhi_};
  
  const double kStepMHT = (l1tmhtemu::kMaxMHT / l1tmhtemu::kMHTBins);
  const double kStepMHTPhi = (l1tmhtemu::kMaxMHTPhi / l1tmhtemu::kMHTPhiBins);

  const unsigned int kPhiBins = 1 << kInternalPhiWidth;

  const float kMaxCosLUTPhi{M_PI};

  const string kLUTdir{"LUTs/"};


  template <typename T>
    T digitizeSignedValue(double value, unsigned int nBits, double lsb) {
    T digitized_value = std::floor(std::abs(value) / lsb);
    T digitized_maximum = (1 << (nBits - 1)) - 1;  // The remove 1 bit from nBits to account for the sign
    if (digitized_value > digitized_maximum)
      digitized_value = digitized_maximum;
    if (value < 0)
      digitized_value = (1 << nBits) - digitized_value;  // two's complement encoding
    return digitized_value;
  }


  //Ouput LUTs to file for use in FW
  template <typename T>
    void writeLUTtoFile(vector<T>(&LUT), const string& filename, const string& delimiter) {
    int fail = system((string("mkdir -p ") + l1tmhtemu::kLUTdir).c_str());
    if (fail)
      throw cms::Exception("BadDir") << __FILE__ << " " << __LINE__ << " could not create directory "
                                     << l1tmhtemu::kLUTdir;

    const string fname = l1tmhtemu::kLUTdir + filename + ".tab";
    ofstream outstream(fname);
    if (outstream.fail())
      throw cms::Exception("BadFile") << __FILE__ << " " << __LINE__ << " could not create file " << fname;

    outstream << "{" << endl;
    for (unsigned int i = 0; i < LUT.size(); i++) {
      if (i != 0)
        outstream << delimiter << endl;
      outstream << LUT[i];
    }
    outstream << endl << "};" << endl;
    outstream.close();
  }


  std::vector<phi_t> generateCosLUT(unsigned int size) {  // Fill cosine LUT with integer values
    float phi = 0;
    std::vector<phi_t> cosLUT;
    for (unsigned int LUT_idx = 0; LUT_idx < size; LUT_idx++) {
      /* cosLUT.push_back((phi_t)(floor(cos(phi) * (kPhiBins - 1)))); */
      cosLUT.push_back(digitizeSignedValue<phi_t>(cos(phi), l1tmhtemu::kInternalPhiWidth, l1tmhtemu::kStepPhi));

      phi += l1tmhtemu::kStepPhi;
    }
    cosLUT.push_back((phi_t)(0));  //Prevent overflow in last bin

    return cosLUT;
  }


  std::vector<MHTphi_t> generateaTanLUT(int cordicSteps) {  // Fill atan LUT with integer values
    std::vector<MHTphi_t> atanLUT;
    for (int cordicStep = 0; cordicStep < cordicSteps; cordicStep++) {
      atanLUT.push_back(MHTphi_t(round((kMHTPhiBins * atan(pow(2,-1*cordicStep))) / (2 * M_PI))));
      /* std::cout<< cordicStep << " -- "  << pow(2,-1*cordicStep)<<"| "; */
    }
    return atanLUT;
  }


  std::vector<Et_t> generatemagNormalisationLUT(int cordicSteps) {
    float val = 1.0;
    std::vector<Et_t> magNormalisationLUT;
    for (int cordicStep = 0; cordicStep < cordicSteps; cordicStep++) {
      val = val / (pow(1 + pow(4, -1*cordicStep), 0.5));
      magNormalisationLUT.push_back(Et_t(round(kMHTBins * val)));
    }
    return magNormalisationLUT;
  }


  struct EtMiss {
    MHT_t Et;
    MHTphi_t Phi;
  };


  EtMiss cordicSqrt(Et_t x, Et_t y, int cordicSteps) {

    std::vector<l1tmhtemu::MHTphi_t> atanLUT;
    atanLUT = generateaTanLUT(cordicSteps);

    std::vector<Et_t> magNormalisationLUT;
    magNormalisationLUT = generatemagNormalisationLUT(cordicSteps);

    Et_t new_x = 0;
    Et_t new_y = 0;

    MHTphi_t phi = 0;
    MHTphi_t new_phi = 0;
    bool sign = false;

    EtMiss ret_etmiss;

    if (x >= 0 && y >= 0) {
      phi = 0;
      sign = true;
      x = x;
      y = y;
    } else if (x < 0 && y >= 0) {
      phi = kMHTPhiBins >> 1;
      sign = false;
      x = -x;
      y = y;
    } else if (x < 0 && y < 0) {
      phi = kMHTPhiBins >> 1;
      sign = true;
      x = -x;
      y = -y;
    } else {
      phi = kMHTPhiBins;
      sign = false;
      x = x;
      y = -y;
    }

    /* std::cout<< " Cordic x: " << x << " Cordic y: " << y << " Cordic phi: " << phi <<endl; */

    for (int step = 0; step < cordicSteps; step++) {
      if (y < 0) {
	new_x = x - (y >> step);
	new_y = y + (x >> step);
      } else {
	new_x = x + (y >> step);
	new_y = y - (x >> step);
      }
    
    if ((y < 0) == sign) {
      new_phi = phi - atanLUT[step];
    } else {
      new_phi = phi + atanLUT[step];
    }

    x = new_x;
    y = new_y;
    phi = new_phi;

    /* std::cout<< " Cordic x: " << x << " Cordic y: " << y << " Cordic phi: " << phi <<std::endl; */

    }

    /* std::cout<< " sqrt (Cordic x * MagNorm / kMHTbins ): " << float(x *  magNormalisationLUT[cordicSteps - 1]) / float(kMHTBins) << " to actual float: " << ( float(x *  magNormalisationLUT[cordicSteps - 1]) / float(kMHTBins) ) * float(kStepPt*kStepPhi) <<std::endl; */

    float sqrtval = ( float(x *  magNormalisationLUT[cordicSteps - 1]) / float(kMHTBins) ) * float(kStepPt*kStepPhi);

    ret_etmiss.Et = std::floor(sqrtval / l1tmhtemu::kStepMHT);
    ret_etmiss.Phi = phi;

    /* std::cout<< " Et: " << ret_etmiss.Et << " phi: " << ret_etmiss.Phi <<std::endl; */

    return ret_etmiss;
  }


}
#endif
