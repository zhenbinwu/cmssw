#ifndef PHASE2GMT_TRACKCONVERTER
#define PHASE2GMT_TRACKCONVERTER

#include "ap_int.h"
#include "L1Trigger/Phase2L1GMT/interface/ConvertedTTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace Phase2L1GMT {

  class TrackConverter {

  public:
    TrackConverter(const edm::ParameterSet& iConfig)  {}
    ~TrackConverter() {}

    std::vector<ConvertedTTTrack> convertTracks(const std::vector<edm::Ptr< l1t::TkMuon::L1TTTrackType > >& tracks) {
      std::vector<ConvertedTTTrack> out;
      for (const auto& t :  tracks)
	out.push_back(convert(t));
      return out;
    }

  private:
    ap_uint<BITSQUALITY> generateQuality(const edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> >& track) {
      return 1;
    }

    ConvertedTTTrack convert(const edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> >& track) {
      ap_uint<1>          charge      = (track->rInv()<0) ? 1 : 0;
      ap_int<BITSTTCURV>  curvature   = track->rInv()*(1<<(BITSTTCURV-1))/maxCurv_;
      ap_int<BITSTTPHI>   phi         = track->phi()*(1<<(BITSTTPHI-1))/maxPhi_;
      ap_int<BITSTTTANL>  tanLambda   = track->tanL()*(1<<(BITSTTTANL-1))/maxTanl_;
      ap_int<BITSTTZ0>    z0          = track->z0()*(1<<(BITSTTZ0-1))/maxZ0_;
      ap_int<BITSTTZ0>    d0          = track->d0()*(1<<(BITSTTD0-1))/maxD0_;
      ap_int<BITSD0>      reducedD0   = ap_int<BITSD0>(d0/2);
      //calculate pt
      ap_uint<BITSTTCURV-1> absCurv  = curvature>0 ? ap_uint<BITSTTCURV-1>(curvature) : ap_uint<BITSTTCURV-1>(-curvature);
      ap_uint<BITSPT> pt = ptLUT[absCurv>>3];
      ap_uint<BITSQUALITY> quality = generateQuality(track);
      ap_uint<BITSTTTANL-1> absTanL = tanLambda >0 ? ap_uint<BITSTTTANL-1>(tanLambda) : ap_uint<BITSTTTANL-1>(-tanLambda);
      ap_uint<BITSETA-1> absEta = etaLUT[absTanL>>4];
      ap_int<BITSETA> eta = tanLambda >0 ?  ap_int<BITSETA>(absEta) : ap_int<BITSETA>(-absEta);

      ConvertedTTTrack convertedTrack(charge,curvature,pt,eta,phi,z0,reducedD0,quality);
      return convertedTrack;
    }
  };
}

#endif



