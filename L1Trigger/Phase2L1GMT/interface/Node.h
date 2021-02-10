#ifndef PHASE2GMT_NODE
#define PHASE2GMT_NODE
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "L1Trigger/Phase2L1GMT/interface/TrackConverter.h"
#include "L1Trigger/Phase2L1GMT/interface/ROITempAssociator.h"
#include "L1Trigger/Phase2L1GMT/interface/TrackMuonMatchAlgorithm.h"

namespace Phase2L1GMT {

  class Node {

  public:
  Node(const edm::ParameterSet& iConfig):
    tt_track_converter_(new TrackConverter(iConfig)),
    roi_assoc_(new ROITempAssociator(iConfig)),
    track_mu_match_(new TrackMuonMatchAlgorithm(iConfig))
  {

  }

  ~Node() 
    {

    }

  std::vector<l1t::TkMuon>  processEvent(const std::vector<edm::Ptr< l1t::TkMuon::L1TTTrackType > >& tracks, 
					 const l1t::ObjectRefBxCollection<l1t::RegionalMuonCand>& muonTracks,
					 const L1TPhase2GMTStubRefVector& stubs) {
    //Convert TT tracks to our internal
    std::vector<ConvertedTTTrack>  convertedTracks =     tt_track_converter_->convertTracks(tracks);

    /* for (const auto& track  : convertedTracks) */
    /*   track.print(); */

    //temporary Associator between stubs and muon tracks till the real stuff arrives from track finders
    std::vector<MuonROI> rois = roi_assoc_->associate(0,muonTracks,stubs);

    
    std::vector<PreTrackMatchedMuon> mu = track_mu_match_->processNonant(convertedTracks,rois);


    std::vector<l1t::TkMuon> out;
    return out;
    }

  private:
    std::unique_ptr<TrackConverter> tt_track_converter_;
    std::unique_ptr<ROITempAssociator> roi_assoc_;
    std::unique_ptr<TrackMuonMatchAlgorithm> track_mu_match_;

  };
}

#endif



