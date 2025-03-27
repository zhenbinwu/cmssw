#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuonEndCapPhase2/interface/EMTFContext.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Utils/DataUtils.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Utils/DebugUtils.h"

#include "L1Trigger/L1TMuonEndCapPhase2/interface/Algo/HybridStubLayer.h"

using namespace emtf::phase2;
using namespace emtf::phase2::algo;

// Members
HybridStubLayer::HybridStubLayer(const EMTFContext& context) : context_(context) {}

// Main function for hybrid stub
void HybridStubLayer::apply (const std::vector<EMTFHitCollection>& bx_window_hits,
                             EMTFTrackCollection& prompt_tracks,
                             EMTFTrackCollection& displaced_tracks) const {
  // ===========================================================================
  // Obtain hybrid stubs variables / IDs
  // ===========================================================================

  // ===========================================================================
  // Matching and form hybrid stubs
  // ===========================================================================
  std::array<EMTFHitCollection, v3::kNumStubsPerSec > hitcollections;

  // Populate
  auto bx_window_hits_rit = bx_window_hits.rbegin();
  auto bx_window_hits_rend = bx_window_hits.rend();

  for (; bx_window_hits_rit != bx_window_hits_rend;
  ++bx_window_hits_rit) {  // Begin loop from latest BX Collection to oldest BX Hit Collection

    const auto& bx_hits = *bx_window_hits_rit;

    for (const auto& hit : bx_hits) {  // Begin loop hits in BX
      //
      // Unpack Hit
      const auto& hit_chamber = hit.emtfChamber();
      const auto& hit_segment = hit.emtfSegment();
      const auto& hit_valid = hit.flagValid();

      emtf_assert(hit_valid);  // segment must be valid

      if( hit_chamber > 113) // Skip ME0 for now
        continue;

      // Calculate Relative BX
      // Note: Uses Hit BX relative to Sector Processor BX
      const auto& hit_bx = hit.bx();
      const int hit_rel_bx = (hit_bx - 0); // Assume only BX=0 now

      // Short-Circuit: Only use Relative BX=0 Segsments
      if (hit_rel_bx != 0) {
        continue;
      }
      hitcollections[hit_chamber % v3::kNumStubsPerSec].push_back(hit);
    }
  }

  hbstub_collection_t stubs;
  for(auto &hits : hitcollections)
  {
    hbstub_collection_t temp = FormHybridPerChamber(hits);
    stubs.insert(stubs.end(), temp.begin(), temp.end());
  }

  // ===========================================================================
  // Matching hybrid stubs to tracks
  // ===========================================================================

  // ===========================================================================
  // Format and sent to GMT
  // ===========================================================================
}


// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::FormHybridPerChamber
//  Description:  
// ===========================================================================
hbstub_collection_t HybridStubLayer::FormHybridPerChamber(EMTFHitCollection& hits ) const
{
  // Sorting by emtfPhi() in ascending order
  std::sort(hits.begin(), hits.end(),
      [](const EMTFHit& a, const EMTFHit& b) { return a.emtfPhi() < b.emtfPhi(); });

   hbstub_collection_t f;
  int segcnt = 0;
  EMTFHit dummy;

  for (unsigned i = 0; i < hits.size(); ++i)
  {
    if (i == hits.size()-1)
    {
      hybrid_stub_t t = CombineHitToHybrid(hits[i], dummy);
      t.stubId = segcnt++;
    };

    int16_t thisphi = hits[i].emtfPhi();
    int16_t nextphi = hits[i+1].emtfPhi();
    if (abs(thisphi - nextphi)< 0.5)
    {
      hybrid_stub_t t = CombineHitToHybrid(hits[i], hits[i+1]);
      t.stubId = segcnt++;
      i++;
    } else
  {
      hybrid_stub_t t = CombineHitToHybrid(hits[i], dummy);
      t.stubId = segcnt++;
    };

  };

  return f;
}       // -----  end of function HybridStubLayer::FormHybridPerChamber  -----


// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::CombineHitToHybrid
//  Description:  /* cursor */
// ===========================================================================
hybrid_stub_t HybridStubLayer::CombineHitToHybrid(EMTFHit& h1, EMTFHit& h2) const
{
  hybrid_stub_t temp;
  temp.valid = h1.flagValid();
  temp.tflayer = lut_.lookup_tflayer(h1.ring(), h1.station());
  temp.etaId = lut_.lookup_etaregion(h1.endcap(), h1.ring());
  temp.secphiId = lut_.lookup_secphi(h1.emtfChamber());

  //Getting Stub 1
  temp.phi_1 = h1.emtfPhi();
  temp.bend_1 = h1.emtfBend();
  temp.theta1_1 = h1.emtfTheta1();
  temp.theta2_1 = h1.emtfTheta2();
  temp.qual1_1 = h1.emtfQual1();
  temp.qual2_1 = h1.emtfQual2();
  temp.time_1 = h1.emtfTime();
  temp.chamber_1 = (h1.emtfChamber()<<2) + h1.emtfSegment();

  if (h2.flagValid())
  {
    //Getting Stub 2
    temp.phi_2 = h2.emtfPhi();
    temp.bend_2 = h2.emtfBend();
    temp.theta1_2 = h2.emtfTheta1();
    temp.theta2_2 = h2.emtfTheta2();
    temp.qual1_2 = h2.emtfQual1();
    temp.qual2_2 = h2.emtfQual2();
    temp.time_2 = h2.emtfTime();
    temp.chamber_2 = (h2.emtfChamber()<<2 )+ h2.emtfSegment();
  }
  else
  {
    temp.phi_2 = 0;
    temp.bend_2 = 0;
    temp.theta1_2 = 0;
    temp.theta2_2 = 0;
    temp.qual1_2 = 0;
    temp.qual2_2 = 0;
  }
  
  return temp;
}       // -----  end of function HybridStubLayer::CombineHitToHybrid  -----


// ===========================================================================
// Obtain hybrid stubs variables / IDs
// ===========================================================================
