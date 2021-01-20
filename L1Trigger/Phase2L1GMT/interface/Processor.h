#ifndef PHASE2GMT_PROCESSOR
#define PHASE2GMT_PROCESSOR

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "L1Trigger/Phase2L1GMT/interface/TrackConverter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


namespace Phase2L1GMT {

  class Processor {

  public:
  Processor(const edm::ParameterSet& iConfig):
    tt_track_converter_(new TrackConverter(iConfig))
  {

  }

  ~Processor() 
    {

    }

    std::vector<l1t::TkMuon>  process(const std::vector<edm::Ptr< l1t::TkMuon::L1TTTrackType > >& tracks) {
      std::vector<ConvertedTTTrack>  convertedTracks =     tt_track_converter_->convertTracks(tracks);
      for (const auto& track  : convertedTracks)
	track.print();
 
      std::vector<l1t::TkMuon> out;
      return out;
    }

  private:
    std::unique_ptr<TrackConverter> tt_track_converter_;
  };
}

#endif



