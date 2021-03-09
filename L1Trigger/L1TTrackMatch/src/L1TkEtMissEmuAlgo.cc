#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuAlgo.h"

namespace L1TkEtMissEmuAlgo {
  std::vector<iglobPhi> FillCosLUT(unsigned int cosLUT_size) {  //Fill cosine LUT with integer values
    float phi = 0;
    std::vector<iglobPhi> cosLUT;
    for (unsigned int LUT_idx = 0; LUT_idx <= cosLUT_size; LUT_idx++) {
      cosLUT.push_back(digitize_Signed<iglobPhi>(cos(phi), 0, 1, N_globPhiBins_));
      phi += (2 * max_TWord_Phi_) / (N_globPhiBins_ - 1);
    }
    return cosLUT;
  }

  // Generate Eta LUT for track to vertex association
  std::vector<iEta> generate_EtaRegions() {
    std::vector<iEta> LUT;
    for (unsigned int q = 0; q < N_etaregions_ + 1; q++) {
      LUT.push_back(digitize_Signed<iEta>(EtaRegions_[q], -max_TWord_Eta_, max_TWord_Eta_, N_etaBins_));
    }
    return LUT;
  }

  std::vector<iZ0> generate_DeltaZBins() {
    std::vector<iZ0> LUT;
    for (unsigned int q = 0; q < N_etaregions_; q++) {
      LUT.push_back(digitize_Signed<iZ0>(
          DeltaZ_[q],
          0,
          max_TWord_Z0_,
          N_z0Bins_ / 2));  //Only half range otherwise deltaZ values midway through integer representation
    }
    return LUT;
  }

}  // namespace L1TkEtMissEmuAlgo