#include <utility>

#include "L1Trigger/L1TMuonEndCapPhase2/interface/Data/HybridLut.h"

using namespace emtf::phase2::data;

// Static
const int HybridLut::kInvalid = -1;

// Member
HybridLut::HybridLut() {
  // TFLayer 3 bits
  // Key: Ring, Station (-1 means take all)
  tflayer_[{1, 1}] = 0;
  tflayer_[{2, 1}] = 4;
  tflayer_[{3, 1}] = 4;
  tflayer_[{1, 2}] = 2;
  tflayer_[{2, 2}] = 2;
  tflayer_[{3, 2}] = 2;
  tflayer_[{1, 3}] = 1;
  tflayer_[{2, 3}] = 1;
  tflayer_[{3, 3}] = 1;
  tflayer_[{1, 4}] = 3;
  tflayer_[{2, 4}] = 3;
  tflayer_[{3, 4}] = 3;

  // Eta region Id 4 bits
  // Key: Ring
  etaRegionId_[1] = 5;
  etaRegionId_[2] = 4;
  etaRegionId_[3] = 3;

  // Phi region ID per sector
  // Key : chamber ID (CSC only)
  // ME11
  phiRegionId_[45] = 0;
  phiRegionId_[0] = 1;
  phiRegionId_[1] = 2;
  phiRegionId_[2] = 3;
  phiRegionId_[9] = 4;
  phiRegionId_[10] = 5;
  phiRegionId_[11] = 6;
  // ME12
  phiRegionId_[46] = 0;
  phiRegionId_[3] = 1;
  phiRegionId_[4] = 2;
  phiRegionId_[5] = 3;
  phiRegionId_[12] = 4;
  phiRegionId_[13] = 5;
  phiRegionId_[14] = 6;
  // ME21
  phiRegionId_[48] = 0;
  phiRegionId_[18] = 1;
  phiRegionId_[19] = 2;
  phiRegionId_[20] = 3;
  // ME22
  phiRegionId_[49] = 0;
  phiRegionId_[21] = 1;
  phiRegionId_[22] = 2;
  phiRegionId_[23] = 3;
  phiRegionId_[24] = 4;
  phiRegionId_[25] = 5;
  phiRegionId_[26] = 6;
  // ME31
  phiRegionId_[50] = 0;
  phiRegionId_[27] = 1;
  phiRegionId_[28] = 2;
  phiRegionId_[29] = 3;
  // ME32
  phiRegionId_[51] = 0;
  phiRegionId_[30] = 1;
  phiRegionId_[31] = 2;
  phiRegionId_[32] = 3;
  phiRegionId_[33] = 4;
  phiRegionId_[34] = 5;
  phiRegionId_[35] = 6;
  // ME41
  phiRegionId_[52] = 0;
  phiRegionId_[36] = 1;
  phiRegionId_[37] = 2;
  phiRegionId_[38] = 3;
  // ME42
  phiRegionId_[53] = 0;
  phiRegionId_[39] = 1;
  phiRegionId_[40] = 2;
  phiRegionId_[41] = 3;
  phiRegionId_[42] = 4;
  phiRegionId_[43] = 5;
  phiRegionId_[44] = 6;
  // ME13
  phiRegionId_[47] = 0;
  phiRegionId_[6] = 1;
  phiRegionId_[7] = 2;
  phiRegionId_[8] = 3;
  phiRegionId_[15] = 4;
  phiRegionId_[16] = 5;
  phiRegionId_[17] = 6;
}

HybridLut::~HybridLut() {
  // Do Nothing
}

void HybridLut::update(const edm::Event&, const edm::EventSetup&) {
 // Do Nothing
}

const int& HybridLut::lookup_tflayer(const int& ring, const int& station) const {
  
  auto found = tflayer_.find(std::make_pair(ring, station));

  if (found == tflayer_.end())
    return HybridLut::kInvalid;

  return found->second;
}

const int HybridLut::lookup_etaregion(const int& endcap, const int& ring) const {
  
  auto found = etaRegionId_.find(ring);

  if (found == etaRegionId_.end())
    return HybridLut::kInvalid;

  return (endcap * found->second);
}

const int& HybridLut::lookup_secphi(const int& chamber) const {
  
  if(chamber > 113)
    return HybridLut::kInvalid;

  auto found = phiRegionId_.find(chamber % 54);

  if (found == phiRegionId_.end())
    return HybridLut::kInvalid;

  return found->second;
}
