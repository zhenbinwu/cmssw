#ifndef L1Trigger_L1TMuonEndCapPhase2_HybridStubLayer_h
#define L1Trigger_L1TMuonEndCapPhase2_HybridStubLayer_h

#include <array>
#include <vector>

#include "L1Trigger/L1TMuonEndCapPhase2/interface/EMTFfwd.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/EMTFTypes.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/EMTFConstants.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Data/HybridLut.h"

namespace emtf::phase2::algo {

  class HybridStubLayer {
    // Static
  private:
    static seg_theta_t calc_theta_median(std::vector<seg_theta_t>);

    // Members
  public:
    HybridStubLayer(const EMTFContext&);

    ~HybridStubLayer() = default;

    void apply (const std::vector<EMTFHitCollection>& bx_window_hits, EMTFTrackCollection& prompt_tracks, EMTFTrackCollection& displaced_tracks) const;

  private:
    const EMTFContext& context_;
    const emtf::phase2::data::HybridLut lut_;

    hybrid_stub_t CombineHitToHybrid(EMTFHit& h1, EMTFHit& h2) const;
    hbstub_collection_t FormHybridPerChamber(EMTFHitCollection& hits ) const;
  };

}  // namespace emtf::phase2::algo

#endif  // L1Trigger_L1TMuonEndCapPhase2_HybridStubLayer_h not defined
