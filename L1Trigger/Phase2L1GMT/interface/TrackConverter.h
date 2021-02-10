#ifndef PHASE2GMT_TRACKCONVERTER
#define PHASE2GMT_TRACKCONVERTER

#include "L1Trigger/Phase2L1GMT/interface/ConvertedTTTrack.h"
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
    uint generateQuality(const edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> >& track) {
      return 1;
    }

    ConvertedTTTrack convert(const edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> >& track) {
      uint                charge      = (track->rInv()<0) ? 1 : 0;
      int                 curvature   = track->rInv()*(1<<(BITSTTCURV-1))/maxCurv_;
      int                 phi         = track->phi()*(1<<(BITSPHI-1))/(M_PI);
      int                 tanLambda   = track->tanL()*(1<<(BITSTTTANL-1))/maxTanl_;
      int                 z0          = track->z0()*(1<<(BITSTTZ0-1))/maxZ0_;
      int                 d0          = track->d0()*(1<<(BITSTTD0-1))/maxD0_;
      int                 reducedD0   = (d0/2);
      //calculate pt
      uint  absCurv  = curvature>0 ? (curvature) : (-curvature);
      uint pt = ptLUT[absCurv>>3];
      uint quality = generateQuality(track);
      uint absTanL = tanLambda >0 ? (tanLambda) : (-tanLambda);
      uint absEta = etaLUT[absTanL>>4];
      int eta = tanLambda >0 ?  (absEta) : (-absEta);

      ConvertedTTTrack convertedTrack(charge,curvature,absEta,pt,eta,phi,z0,reducedD0,quality);
      convertedTrack.setOfflineQuantities(track->momentum().transverse(),track->eta(),track->phi());
      return convertedTrack;
    }
  };
}

#endif



