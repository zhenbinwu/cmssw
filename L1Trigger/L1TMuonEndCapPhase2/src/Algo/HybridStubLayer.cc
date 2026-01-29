#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "L1Trigger/L1TMuonEndCapPhase2/interface/EMTFContext.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Utils/DataUtils.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Utils/DebugUtils.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Utils/TPUtils.h"

#include "L1Trigger/L1TMuonEndCapPhase2/interface/Algo/HybridStubLayer.h"

using namespace emtf::phase2;
using namespace emtf::phase2::algo;

// Members
HybridStubLayer::HybridStubLayer(const EMTFContext& context) : context_(context) {
  emtfhitFile = new std::ofstream("emtfhit.txt");
  hybridFile = new std::ofstream("hybrid.txt");
}

// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::ConvertPhiBits
//  Description:  
// ===========================================================================
int HybridStubLayer::ConvertPhiBits(int phi) const
{
  const int LOCPhiBits = 10;
  const float LSBPhi = 2 *M_PI / pow(2, LOCPhiBits);

  float phirad = emtf::phase2::tp::calcPhiLocRadFromInt(phi);
  int phiBit = round(phirad  /  LSBPhi);

  return phiBit;
}       // -----  end of function HybridStubLayer::ConvertPhiBits  -----

// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::ConvertPhiBits
//  Description:  
// ===========================================================================
int HybridStubLayer::ConvertPhiBits(int sector, int phi) const
{
  const int LOCPhiBits = 10;
  const float LSBPhi = 2 *M_PI / pow(2, LOCPhiBits);

  float phirad = emtf::phase2::tp::calcPhiGlobRadFromLoc(sector, emtf::phase2::tp::calcPhiLocRadFromInt(phi));
  int phiBit = round(phirad  /  LSBPhi);
  std::cout << " ================ sector " << sector << " phi " << phi <<" global phi "<< phirad <<" phiBit " << phiBit << std::endl;

  return phiBit;
}       // -----  end of function HybridStubLayer::ConvertPhiBits  -----

// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::ConvertEtaBits
//  Description:  
// ===========================================================================
int HybridStubLayer::ConvertEtaBits(int endcap, int theta) const
{
  const int LOCEtaBits = 8;
  const float LSBeta = 2 *M_PI / pow(2, LOCEtaBits);

  float track_theta = emtf::phase2::tp::calcThetaRadFromInt(theta);
  float track_eta = -1 * std::log(std::tan(track_theta / 2));
  track_eta *= endcap;

  int etabit = round(track_eta / LSBeta);
  std::cout << " global theta " << theta <<" rad theta " << track_theta << " eta? " << track_eta << " etabit " << etabit << std::endl;

  return etabit;
}       // -----  end of function HybridStubLayer::ConvertEtaBits  -----

// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::ConvertEta
//  Description:  
// ===========================================================================
float HybridStubLayer::ConvertEta(int endcap, float theta) const
{
  float track_theta = emtf::phase2::tp::degToRad(theta);
  float track_eta = -1 * std::log(std::tan(track_theta / 2));
  //track_eta *= float(endcap);
  return track_eta;
}       // -----  end of function HybridStubLayer::ConvertEtaBits  -----

// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::PrintEMTFHit
//  Description:  
// ===========================================================================
bool HybridStubLayer::PrintEMTFHit(const EMTFHit &hit) const
{
  if (emtfhitFile-> is_open())
  {
    *emtfhitFile <<  hit.endcap() <<" "
      << hit.sector() <<" " 
      << hit.ring() <<" " 
      << hit.station() <<" " 
      << hit.emtfChamber() <<" " 
      << hit.emtfSegment() <<" " 
      << hit.emtfPhi() <<" " 
      << hit.emtfTheta1() <<" " 
      << hit.emtfTheta2() <<" " 
      << std::endl;
  }

  return true;
}       // -----  end of function HybridStubLayer::PrintEMTFHit  -----

// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::PrintHybridStub
//  Description:  
// ===========================================================================
bool HybridStubLayer::PrintHybridStub(const l1t::MuonStub& s ) const
{
  if (hybridFile->is_open())
  {
    *hybridFile << s.etaRegion() <<" "
      << s.phiRegion() <<" "
      << s.depthRegion() << " "
      << s.tfLayer() << " "
      << s.id() << " "
      << s.coord1()<<" "
      << s.coord2()<<" "
      << s.eta1()<<" "
      << s.eta2()<<" " << s.quality() << " "
      << s.etaQuality()<< " " << s.type() << std::endl;
  }
  return true;
}       // -----  end of function HybridStubLayer::PrintHybridStub  -----

// Main function for hybrid stub
void HybridStubLayer::apply ( const std::vector<EMTFHitCollection>& bx_window_hits,
    EMTFTrackCollection& prompt_tracks,
    EMTFTrackCollection& displaced_tracks, 
    l1t::MuonStubCollection& out_stubs, 
    l1t::SAMuonCollection& out_samuons) const {
  // ===========================================================================
  // Obtain hybrid stubs variables / IDs
  // ===========================================================================

  // ===========================================================================
  // Matching and form hybrid stubs
  // ===========================================================================
  for(auto &v : hitcollections)
  {
    v.clear();
  }

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
      //const auto& hit_segment = hit.emtfSegment();
      const auto& hit_valid = hit.flagValid();

      emtf_assert(hit_valid);  // segment must be valid

      // TODO: Skip ME0 for now
      if( hit_chamber > 113) 
        continue;

      // TODO: Skip GE11 for too?
      if(( hit_chamber == 99) || ( hit_chamber >= 54 && hit_chamber <= 56 ) 
          ||  ( hit_chamber >= 63 && hit_chamber <= 65 ) )
        continue;

      // Skip the overlap region
      if(( hit_chamber >= 45 && hit_chamber <= 53 ) ||  ( hit_chamber >= 99 && hit_chamber <= 108 ) )
        continue;

      PrintEMTFHit(hit);

      //// TODO: Skip ME11 in high thera (belwo strip cut) 
      //if( hit.ring() > 3) 
        //continue;

      // Calculate Relative BX
      // Note: Uses Hit BX relative to Sector Processor BX
      const auto& hit_bx = hit.bx();
      const int hit_rel_bx = (hit_bx - 0); // Assume only BX=0 now

      // Short-Circuit: Only use Relative BX=0 Segsments
      if (hit_rel_bx != 0) {
        continue;
      }

      // Put hit into collection per chamber ID
      //hitcollections[2].push_back(hit);
      hitcollections[hit_chamber % v3::kNumStubsPerSec].push_back(hit);
    }
  }

  hbstub_collection_t stubs;
  for (int i = 0; i < v3::kNumStubsPerSec; ++i)
  {
    //hbstub_collection_t temp;
    hbstub_collection_t temp = FormHybridPerChamber(i, hitcollections.at(i));
    stubs.insert(stubs.end(), temp.begin(), temp.end());
  }

  // ===========================================================================
  // Matching hybrid stubs to tracks
  // ===========================================================================
  //MatchHybridStubsToTracks(stubs, prompt_tracks, displaced_tracks);

  // ===========================================================================
  // Format and sent to GMT
  // ===========================================================================
  for(auto &s : stubs)
  {
    //int etaRegion = 0;
    //int phiRegion = 0;
    //uint tflayer = 0;
    int qual = 1;
    if (s.phi_2 != 0)
    {
      qual |= 0x1;
    }

    l1t::MuonStub t(s.etaId, s.globphiId, s.depthId, s.tflayer, 
        int(s.phi_1), int(s.phi_2), s.chamber_1, 0, qual,  
        int(s.theta1_1), int(s.theta1_2), 1, 0);

    t.setOfflineQuantities(s.glob_phi_1, s.glob_phi_2, s.glob_theta_1, s.glob_theta_2);

    out_stubs.push_back(t);
    PrintHybridStub(t);
  }
}

// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::FormHybridPerChamber
//  Description:  Forms hybrid stubs for a given chamber from its hits.
//                This function is called per chamber by HybridStubLayer::apply.
// ===========================================================================
hbstub_collection_t HybridStubLayer::FormHybridPerChamber(const int chamberID, const EMTFHitCollection& chamber_hits) const {

  hbstub_collection_t chamber_stubs;

  // TODO: Implement the logic to form hybrid stubs from the input 'chamber_hits'.
  // This function should process the hits within a single chamber and
  // return a collection of 'hybrid_stub_t' objects.
  if (chamberID > 0 ) // dummy for now
  {
    chamber_stubs = BuildCSCRPCStubs(chamber_hits);
  }

  // For now, returning an empty collection.
  return chamber_stubs;
}

// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::MatchHybridStubsToPromptTrack
//  Description:  Matches hybrid stubs to a single prompt EMTF track
// ===========================================================================
void HybridStubLayer::MatchHybridStubsToPromptTrack(const hybrid_stub_t& stub, 
                                                    EMTFTrack& track,
                                                    const hbstub_collection_t& stubs) const {
  /*
  // Skip low-quality stubs
  //if (stub.quality() < kMinStubQuality) {
    //return;
  //}

  const auto& stub_bx = stub.bx();
  const auto stub_hit_id = getStubHitId(stub);

  // Validate hit ID
  if (!isValidHitId(stub_hit_id)) {
    return;
  }

  // Check BX compatibility
  if (track.bx() != stub_bx) {
    return;
  }

  // Check if stub's hit ID matches any of the track's site hits
  bool match_found = false;
  for (int i = 0; i < kMaxSiteHits; ++i) {
    if (track.site_mask()[i] != 0) {  // Only check if site is active
      if (track.site_hits()[i] == stub_hit_id) {
        match_found = true;
        break;
      }
    }
  }

  // If match found and track has space, add stub
  if (match_found && track.stubs().size() < kMaxStubsPerTrack) {
    track.addStub(stub);
    
    // Log successful match
    edm::LogVerbatim("L1TEMTFpp")
      << "Matched stub to track:"
      << " track_id=" << track.id()
      << " stub_id=" << stub_hit_id
      << " phi=" << stub.emtfPhi()
      << " quality=" << stub.quality();
  }  
  */
}

// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::MatchHybridStubsToTracks
//  Description:  Matches hybrid stubs to EMTF tracks based on phi angle and other parameters
// ===========================================================================
void HybridStubLayer::MatchHybridStubsToTracks(const hbstub_collection_t& stubs, 
                                               EMTFTrackCollection& prompt_tracks, 
                                               EMTFTrackCollection& displaced_tracks) const {

  /* 
  // Helper function to get stub hit ID
  auto getStubHitId = [](const hybrid_stub_t& stub) -> uint16_t {
    // Extract hit ID from stub - this is implementation-specific
    // For now, we'll assume it's stored in the stub's data structure
    return stub.emtfChamber();  // This is just a placeholder
  };

  // Helper function to validate hit ID
  bool isValidHitId(uint16_t hit_id) const {
    if (hit_id == 0 || hit_id >= v3::kNumChambers) {
      edm::LogWarning("L1TEMTFpp") << "Invalid hit ID: " << hit_id;
      return false;
    }
    return true;
  }

  // Process prompt tracks first
  for (auto& track : prompt_tracks) {
    // Skip low-quality tracks
    if (track.quality() < kMinTrackQuality) {
      continue;
    }

    // Process each stub
    for (const auto& stub : stubs) {
      MatchHybridStubsToPromptTrack(stub, track, stubs);
    }
  }

  // Process displaced tracks
  for (auto& track : displaced_tracks) {
    // Skip low-quality tracks
    if (track.quality() < kMinTrackQuality) {
      continue;
    }

    // Process each stub
    for (const auto& stub : stubs) {
      MatchHybridStubsToPromptTrack(stub, track, stubs);
    }
  }
      // Skip low-quality stubs
      if (stub.quality() < kMinStubQuality) {
        continue;
      }

      const auto& stub_bx = stub.bx();
      const auto stub_hit_id = getStubHitId(stub);

      // Validate hit ID
      if (!isValidHitId(stub_hit_id)) {
        continue;
      }

      // Check BX compatibility
      if (track_bx != stub_bx) {
        continue;
      }

      // Check if stub's hit ID matches any of the track's site hits
      bool match_found = false;
      for (int i = 0; i < kMaxSiteHits; ++i) {
        if (track_site_mask[i] != 0) {  // Only check if site is active
          if (track_site_hits[i] == stub_hit_id) {
            match_found = true;
            break;
          }
        }
      }

      // If match found and track has space, add stub
      if (match_found && track.stubs().size() < kMaxStubsPerTrack) {
        track.addStub(stub);
        
        // Log successful match
        edm::LogVerbatim("L1TEMTFpp")
          << "Matched stub to track:"
          << " track_id=" << track.id()
          << " stub_id=" << stub_hit_id
          << " phi=" << stub.emtfPhi()
          << " quality=" << stub.quality();
      }
    }
  }

  // Process displaced tracks
  for (auto& track : displaced_tracks) {
      const auto& stub_bx = stub.bx();
      const auto stub_hit_id = getStubHitId(stub);

      // Check BX compatibility
      if (track_bx != stub_bx) continue;

      // Check if stub's hit ID matches any of the track's site hits
      bool match_found = false;
      for (int i = 0; i < kMaxSiteHits; ++i) {
        if (track_site_mask[i] != 0) {  // Only check if site is active
          if (track_site_hits[i] == stub_hit_id) {
            match_found = true;
            break;
          }
        }
      }

      // If match found and track has space, add stub
      if (match_found && track.stubs().size() < kMaxStubsPerTrack) {
        track.addStub(stub);
      }
    }
  }

  // Process displaced tracks
  for (auto& track : displaced_tracks) {
    const auto& track_bx = track.bx();
    const auto& track_site_hits = track.site_hits();
    const auto& track_site_mask = track.site_mask();

    for (const auto& stub : stubs) {
      const auto& stub_phi = stub.emtfPhi();
      const auto& stub_bx = stub.bx();

      if (track_bx != stub_bx) continue;
      if (std::abs(track_phi - stub_phi) > kPhiWindow) continue;

        if (track.stubs().size() < kMaxStubsPerTrack) {
        track.addStub(stub);
      }
    }
  }
*/
}

// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::BuildCSCRPCStubs
//  Description:  Build stubs from CSC and RPC hits
// ===========================================================================
// Helper function to find the best match for a hit within the phi cut window
std::pair<int, float> HybridStubLayer::FindBestMatch(const EMTFHit& hit, const EMTFHitCollection& hits, std::vector<bool>& used) const
{
  int16_t thisphi = ConvertPhiBits(hit.emtfPhi());
  unsigned best_match = -1;
  float best_dphi = 999999;
  std::vector<std::pair<int, float> > matches;

  // Find the best match within dphi_cut
  for (unsigned j = 0; j < hits.size(); ++j) {
    if (used[j] ) continue;
    if ( (hits[j].emtfPhi() == hit.emtfPhi()) &&
        (hits[j].emtfTheta1() == hit.emtfTheta1()) &&
        (hits[j].emtfTheta2() == hit.emtfTheta2()) ) continue;

    // within same sector, no need to wrap
    float dphi = std::abs(thisphi - ConvertPhiBits(hits[j].emtfPhi()));
    const EMTFHit& h = hits[j];
    std::cout <<" compared h endcap " << h.endcap()<<", sector  "<< h.sector()<<", chamber " << h.chamber() 
      <<",  CSC chamber" << h.chamber() % v3::kNumStubsPerSec << " , type ?" << h.subsystem() 
      <<",  phi " << h.emtfPhi() << std::endl;
    if (dphi < dphi_cut_)
    {
      matches.push_back(std::make_pair(j, dphi));
    }
    if (dphi == 0)
    {
      std::cout << " not the same hit but dphi ==0? " << std::endl;
    }
    //std::cout << " this phi " << thisphi <<" newphi" << hits[j].emtfPhi() << " deltaphi "<< dphi <<" cut "<< dphi_cut_ << std::endl;
    if (dphi < dphi_cut_ && dphi < best_dphi) {
      best_dphi = dphi;
      best_match = j;
    }
  }

  if (matches.size()>1)
  {
		for(auto [i, j] : matches)
		{
			std::cout << "++++++ in matches " << i <<" ; " << j  << std::endl;
		}
    // Find the pair with the maximum second element
    auto max_second_it = std::max_element(matches.begin(), matches.end(),
                                          [](const auto& p1, const auto& p2) {
                                              return p1.second < p2.second;
                                          });

    if (max_second_it != matches.end()) {
        std::cout << "Max pair (by second element): (" << max_second_it->first << ", " << max_second_it->second << ")" << std::endl;
				best_match = max_second_it->first;
				best_dphi = max_second_it->second;

				for(auto [i, j] : matches)
				{
					used[i]=true;
				}
    } else {
        std::cout << "Vector is empty." << std::endl;
    }
  }
  return {best_match, best_dphi};
}

hbstub_collection_t HybridStubLayer::BuildCSCRPCStubs(const EMTFHitCollection& hits_ ) const
{
  EMTFHitCollection hits = hits_;
  // Sort hits by phi angle
  std::sort(hits.begin(), hits.end(),
      [](EMTFHit a, EMTFHit b) { return a.emtfPhi() <= b.emtfPhi(); });

  hbstub_collection_t f;
  int segcnt = 0;
  EMTFHit dummy;
  std::vector<bool> used(hits.size(), false);  // Track which hits are used

  // Process each hit
  for (unsigned i = 0; i < hits.size(); ++i) {
    if (used[i]) continue;  // Skip already used hits

    const auto& hit = hits[i];

    auto [best_match, best_dphi] = FindBestMatch(hit, hits, used);

    std::cout <<" target hit endcap " << hit.endcap()<<", sector  "<< hit.sector()<<", chamber " << hit.chamber() 
      <<",  CSC chamber" << hit.chamber() % v3::kNumStubsPerSec << " , type ?" << hit.subsystem() 
      <<",  phi " << hit.emtfPhi() << ", best match " << best_match <<", best dphi " << best_dphi  << std::endl;

    if (best_dphi == 0)
    {
      std::cout << "--------+++++++++++++++++++++++" << std::endl;
      std::cout << "this hit " << hit.emtfChamber() <<" ; " << hit.emtfSegment() <<" ; " 
        << hit.emtfTheta1() << " ; " << hit.emtfTheta2()<< std::endl;
      std::cout << "that hit " << hits[best_match].emtfChamber() <<" ; " << hits[best_match].emtfSegment() 
        <<" ; " << hits[best_match].emtfTheta1() << " ; " << hits[best_match].emtfTheta2()<< std::endl;
    }

    hybrid_stub_t t;

    // Combine hits if we found a good match
    if (best_match != -1) {
      assert(best_match >=0 && best_match < static_cast<int>( hits.size()));
      t = CombineHitToHybrid(hit, hits[best_match]);
      t.stubId = segcnt++;
      used[i] = true;
      used[best_match] = true;
    } else {
      t = CombineHitToHybrid(hit, dummy);
      t.stubId = segcnt++;
      used[i] = true;
    }
    f.push_back(t);
  }

  return f;
}       // -----  end of function HybridStubLayer::BuildCSCRPCStubs  -----


// ===  FUNCTION  ============================================================
//         Name:  HybridStubLayer::CombineHitToHybrid
//  Description:  /* cursor */
// ===========================================================================
hybrid_stub_t HybridStubLayer::CombineHitToHybrid(const EMTFHit& h1,const  EMTFHit& h2) const
{

  hybrid_stub_t temp;
  temp.valid = h1.flagValid();
  int ring = h1.ring();
  if (ring == 4)
  {
    ring = 1;
  }
  temp.tflayer = lut_.lookup_tflayer(ring, h1.station());
  temp.etaId = lut_.lookup_etaregion(h1.endcap(), h1.ring());
  temp.secphiId = lut_.lookup_secphi(h1.emtfChamber());
  //temp.secphiId = h1.chamber();
  temp.globphiId = lut_.lookup_phi(h1.sector(), h1.emtfChamber());
  temp.depthId = h1.station();

  //Getting Stub 1
  temp.phi_1 = ConvertPhiBits(h1.sector(), h1.emtfPhi());
  temp.bend_1 = h1.emtfBend();
  temp.theta1_1 = ConvertEtaBits(h1.endcap(), h1.emtfTheta1());
  temp.theta2_1 = ConvertEtaBits(h1.endcap(), h1.emtfTheta2());
  temp.qual1_1 = h1.emtfQual1();
  temp.qual2_1 = h1.emtfQual2();
  temp.time_1 = h1.emtfTime();
  temp.chamber_1 = (h1.emtfChamber()<<2) + h1.emtfSegment();
  temp.glob_phi_1 = emtf::phase2::tp::degToRad(h1.globPhi());
  temp.glob_theta_1 = ConvertEta(h1.endcap(), h1.globTheta());
  temp.glob_perp_1 = h1.globPerp();
  temp.glob_z_1 = h1.globZ();
  temp.glob_time_1 = h1.globTime();

  if (h2.flagValid())
  {
    //Getting Stub 2
    temp.phi_2 = ConvertPhiBits(h2.sector(), h2.emtfPhi());
    temp.bend_2 = h2.emtfBend();
    temp.theta1_2 = ConvertEtaBits(h2.endcap(), h2.emtfTheta1());
    temp.theta2_2 = ConvertEtaBits(h2.endcap(), h2.emtfTheta2());
    temp.qual1_2 = h2.emtfQual1();
    temp.qual2_2 = h2.emtfQual2();
    temp.time_2 = h2.emtfTime();
    temp.chamber_2 = (h2.emtfChamber()<<2 )+ h2.emtfSegment();
    temp.glob_phi_2 = emtf::phase2::tp::degToRad(h1.globPhi());
    temp.glob_theta_2 = ConvertEta(h2.endcap(), h2.globTheta());
    temp.glob_perp_2 = h2.globPerp();
    temp.glob_z_2 = h2.globZ();
    temp.glob_time_2 = h2.globTime();
    temp.qual &= 0x1;
  } else {
    temp.phi_2 = 0;
    temp.bend_2 = 0;
    temp.theta1_2 = 0;
    temp.theta2_2 = 0;
    temp.qual1_2 = 0;
    temp.qual2_2 = 0;
    temp.glob_phi_2 = 0;
    temp.glob_theta_2 = 0;
    temp.glob_perp_2 = 0;
    temp.glob_z_2 = 0;
    temp.glob_time_2 = 0;
  }
  
  return temp;
}       // -----  end of function HybridStubLayer::CombineHitToHybrid  -----


// ===========================================================================
// Obtain hybrid stubs variables / IDs
// ===========================================================================
