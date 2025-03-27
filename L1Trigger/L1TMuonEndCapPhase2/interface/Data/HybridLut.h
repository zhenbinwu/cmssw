#ifndef L1Trigger_L1TMuonEndCapPhase2_HybridLut_h
#define L1Trigger_L1TMuonEndCapPhase2_HybridLut_h

#include <map>
#include <tuple>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

namespace emtf::phase2::data {

  class HybridLut {
    // Static
  public:
    static const int kInvalid;

    // Member
  public:
    HybridLut();

    ~HybridLut();

    void update(const edm::Event&, const edm::EventSetup&);

  const int& lookup_tflayer(const int& ring, const int& station) const;
  const int lookup_etaregion(const int& endcap, const int& ring) const;
  const int& lookup_secphi(const int& chamber) const;

  private:
    // Key: Subsystem, Station, Ring
    // Value: Hybrid
    std::map<std::tuple<int, int, int>, int> lut_;

    // TFLayer 3 bits
    // Key: Ring, Station
    std::map<std::pair<int, int>, int> tflayer_;

    // Eta region ID 4 bits
    // Key: Ring
    std::map<int, int> etaRegionId_;

    // Phi region ID 6 bits for 360 degrees, so 10 value per chamber
    // Key: Ring
    std::map<int, int> phiRegionId_;
  };

}  // namespace emtf::phase2::data

#endif  // L1Trigger_L1TMuonEndCapPhase2_HybridLut_h
