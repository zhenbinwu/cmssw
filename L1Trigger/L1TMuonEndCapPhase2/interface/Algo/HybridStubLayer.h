#ifndef L1Trigger_L1TMuonEndCapPhase2_HybridStubLayer_h
#define L1Trigger_L1TMuonEndCapPhase2_HybridStubLayer_h

#include <array>
#include <vector>
#include <fstream>

#include "L1Trigger/L1TMuonEndCapPhase2/interface/EMTFfwd.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/EMTFTypes.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/EMTFConstants.h"
#include "L1Trigger/L1TMuonEndCapPhase2/interface/Data/HybridLut.h"
#include "DataFormats/L1TMuonPhase2/interface/MuonStub.h"
#include "DataFormats/L1TMuonPhase2/interface/SAMuon.h"

namespace emtf::phase2::algo {

  class HybridStubLayer {
    // Static
  private:
    static seg_theta_t calc_theta_median(std::vector<seg_theta_t>);
    // Constants for matching
    const int kMaxStubsPerTrack = 3;  // Maximum number of stubs per track
    const int kMaxSiteHits = 12;  // Maximum number of site hits to check
    const int kMinTrackQuality = 3;  // Minimum track quality for matching
    const int kMinStubQuality = 2;   // Minimum stub quality for matching

    std::ofstream *emtfhitFile;
    std::ofstream *hybridFile;

    // Members
  public:
    HybridStubLayer(const EMTFContext&);

    ~HybridStubLayer() = default;

    void apply ( const std::vector<EMTFHitCollection>& bx_window_hits, EMTFTrackCollection& prompt_tracks, EMTFTrackCollection& displaced_tracks, l1t::MuonStubCollection& out_stubs, l1t::SAMuonCollection& out_samuons) const;

  private:
    const EMTFContext& context_;
    const emtf::phase2::data::HybridLut lut_;
    mutable std::vector<bool> used;
    mutable std::array<EMTFHitCollection, v3::kNumStubsPerSec > hitcollections;
    float dphi_cut_ = 9;  // Maximum phi difference for combining hits into a stub
    // 10*LSB = ~0.06rad
    // 9*LSB = ~0.05rad

    // Find the best match for a hit within the phi cut window
    std::pair<int, float> FindBestMatch(const EMTFHit& hit, const EMTFHitCollection& hits, std::vector<bool>& used) const;

    hbstub_collection_t FormHybridPerChamber(const int chamberID, const EMTFHitCollection& chamber_hits) const;

    hybrid_stub_t CombineHitToHybrid(const EMTFHit& h1,const EMTFHit& h2) const;
    hbstub_collection_t BuildCSCRPCStubs(const EMTFHitCollection& hits ) const;
    void MatchHybridStubsToTracks(const hbstub_collection_t& stubs, EMTFTrackCollection& prompt_tracks, EMTFTrackCollection& displaced_tracks) const;

    // Helper function for matching individual stubs to tracks
    void MatchHybridStubsToPromptTrack(const hybrid_stub_t& stub, 
                                       EMTFTrack& track,
                                       const hbstub_collection_t& stubs) const;
    int ConvertPhiBits(int phi) const;
    bool PrintEMTFHit(const EMTFHit &hit) const;
    bool PrintHybridStub(const l1t::MuonStub& s ) const;
    int ConvertPhiBits(int sector, int phi) const;
    int ConvertEtaBits(int endcap, int theta) const;
    float ConvertEta(int endcap, float theta) const;
  };

}  // namespace emtf::phase2::algo

#endif  // L1Trigger_L1TMuonEndCapPhase2_HybridStubLayer_h not defined
