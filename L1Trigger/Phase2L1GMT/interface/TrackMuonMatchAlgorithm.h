#ifndef PHASE2GMT_TRACKMUONMATCHALGO
#define PHASE2GMT_TRACKMUONMATCHALGO

#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TMuonPhase2/interface/MuonStub.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1Trigger/interface/L1TObjComparison.h"
#include "L1Trigger/Phase2L1GMT/interface/TrackConverter.h"
#include "L1Trigger/Phase2L1GMT/interface/PreTrackMatchedMuon.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "L1Trigger/Phase2L1GMT/interface/Constants.h"

namespace Phase2L1GMT {
  
  const int BITSPROP_C1 = 10;
  const int BITSPROP_C2 = 10;
  const int BITSPROP_SIGMAC1_A = 10;
  const int BITSPROP_SIGMAC1_B = 10;
  const int BITSPROP_SIGMAC2_A = 10;
  const int BITSPROP_SIGMAC2_B = 10;
  const int BITSPROP_SIGMAETA1_A = 10;
  const int BITSPROP_SIGMAETA_B = 10;
  const int BITSPROP_SIGMAETA2_A = 10;




  typedef struct {
    ap_int<BITSSTUBCOORD1> coord1;
    ap_uint<BITSSTUBCOORD1-1> sigma_coord1;
    ap_int<BITSSTUBCOORD2> coord2;
    ap_uint<BITSSTUBCOORD2-1> sigma_coord2;
    ap_int<BITSSTUBETA> eta;
    ap_uint<BITSSTUBETA-1> sigma_eta1;
    ap_uint<BITSSTUBETA-1> sigma_eta2;
    ap_uint<1> valid;
    ap_uint<1> is_barrel;
  } propagation_t;

  typedef struct {
    ap_uint<BITSMATCHQUALITY> quality;
    ap_uint<BITSSTUBID> id;
    ap_uint<1> valid;
    l1t::RegionalMuonCandRef muRef;
    l1t::MuonStubRef stubRef;
  } match_t;

  class TrackMuonMatchAlgorithm {

  public:
  TrackMuonMatchAlgorithm(const edm::ParameterSet& iConfig):
    verbose_(1)
    {

    }

  ~TrackMuonMatchAlgorithm() 
    {

    }


  std::vector<PreTrackMatchedMuon>  processNonant(const std::vector<ConvertedTTTrack>& convertedTracks,
						  const std::vector<MuonROI>& rois) {

    std::vector<PreTrackMatchedMuon> preMuons;
    for (const auto& track : convertedTracks) {
      PreTrackMatchedMuon mu = processTrack(track,rois);
      if (mu.quality()>0 && preMuons.size()<32)
	preMuons.push_back(mu);
    }
    std::vector<PreTrackMatchedMuon> cleanedMuons  = clean(preMuons);
    return cleanedMuons;

  }

  std::vector<PreTrackMatchedMuon> cleanNeighbor(const std::vector<PreTrackMatchedMuon>&  muons,const std::vector<PreTrackMatchedMuon>&  muonsPrevious,const std::vector<PreTrackMatchedMuon>&  muonsNext, bool equality) {
    std::vector<PreTrackMatchedMuon> out;
    if (muons.size()==0)
      return out;

    if (verbose_==1) {
      printf("-----Cleaning Up Muons in the neighbours\n");
      printf("Before:\n");
    }


    for (uint i=0;i<muons.size();++i) {
      if (verbose_==1)
	muons[i].print();

      bool keep=true;
      for (uint j=0;j<muonsPrevious.size();++j) {
	bool overlap=false;
	if (muons[i].stubID0()==muons[j].stubID0() && muons[i].stubID0()!=511)
	  overlap=true;
	if (muons[i].stubID1()==muons[j].stubID1() && muons[i].stubID1()!=511)
	  overlap=true;
	if (muons[i].stubID2()==muons[j].stubID2() && muons[i].stubID2()!=511)
	  overlap=true;
	if (muons[i].stubID3()==muons[j].stubID3() && muons[i].stubID3()!=511)
	  overlap=true;
	if (muons[i].stubID4()==muons[j].stubID4() && muons[i].stubID4()!=511)
	  overlap=true;

	if ((overlap && muons[i].quality()<muons[j].quality() && !equality)||(overlap && muons[i].quality()<=muons[j].quality() &&equality))
	  keep=false;
      }
      for (uint j=0;j<muonsNext.size();++j) {
	bool overlap=false;
	if (muons[i].stubID0()==muons[j].stubID0() && muons[i].stubID0()!=511)
	  overlap=true;
	if (muons[i].stubID1()==muons[j].stubID1() && muons[i].stubID1()!=511)
	  overlap=true;
	if (muons[i].stubID2()==muons[j].stubID2() && muons[i].stubID2()!=511)
	  overlap=true;
	if (muons[i].stubID3()==muons[j].stubID3() && muons[i].stubID3()!=511)
	  overlap=true;
	if (muons[i].stubID4()==muons[j].stubID4() && muons[i].stubID4()!=511)
	  overlap=true;

	if ((overlap && muons[i].quality()<muons[j].quality() && !equality)||(overlap && muons[i].quality()<=muons[j].quality() &&equality))
	  keep=false;
      }



      if (keep) {
	if (verbose_==1)
	  printf("kept\n");
	out.push_back(muons[i]);
      }
      else {
	if (verbose_==1)
	  printf("discarded\n");


      }
    }
    return out;
  }


  std::vector<l1t::TrackerMuon> convert(std::vector<PreTrackMatchedMuon>&  muons,uint maximum) {
    std::vector<l1t::TrackerMuon> out;
    for (const auto& mu : muons)  {
      if (out.size()==maximum)
	break;
      l1t::TrackerMuon muon(mu.trkPtr(),mu.charge(),mu.pt(),mu.eta(),mu.phi(),mu.z0(),mu.d0(),mu.quality());			    
      muon.setMuonRef(mu.muonRef());
      for (const auto& stub: mu.stubs())
	muon.addStub(stub);
      out.push_back(muon);
      if(verbose_==1) {
	printf("Final Muon:");
	muon.print();
      }
    }
    return out;
  }    


  private:

  int verbose_;

  propagation_t propagate(const ConvertedTTTrack& track, uint layer) {
    ap_uint<BITSPROP_C1>          prop_coord1 = ap_uint<BITSPROP_C1>(0);
    ap_uint<BITSPROP_C2>          prop_coord2 = 0;
    ap_uint<BITSPROP_SIGMAC1_A>   res0_coord1 = 0;
    ap_uint<BITSPROP_SIGMAC1_B>   res1_coord1= 0;
    ap_uint<BITSPROP_SIGMAC2_A>   res0_coord2 = 0;
    ap_uint<BITSPROP_SIGMAC2_B>   res1_coord2 = 0;
    ap_uint<BITSPROP_SIGMAETA1_A> res0_eta1 = 0;
    ap_uint<BITSPROP_SIGMAETA_B>  res1_eta = 0;
    ap_uint<BITSPROP_SIGMAETA2_A> res0_eta2 = 0;
    ap_uint<1>                    valid = ap_uint<1>(0);
    ap_uint<1>                    is_barrel = 0;

    uint reducedAbsEta = track.abseta()/32;

    if (layer ==0 ) {
      prop_coord1        =      lt_prop_coord1_0[reducedAbsEta];
      prop_coord2        =      lt_prop_coord2_0[reducedAbsEta];
      res0_coord1        =      lt_res0_coord1_0[reducedAbsEta];
      res1_coord1        =      lt_res1_coord1_0[reducedAbsEta];
      res0_coord2        =      lt_res0_coord2_0[reducedAbsEta];
      res1_coord2        =      lt_res1_coord2_0[reducedAbsEta];
      res0_eta1          =      lt_res0_eta1_0[reducedAbsEta];
      res1_eta           =      lt_res1_eta_0[reducedAbsEta];
      res0_eta2          =      lt_res0_eta2_0[reducedAbsEta];
      valid              =      lt_prop_coord1_0_valid[reducedAbsEta];
      is_barrel          =      reducedAbsEta<barrelLimit0_ ? 1:0;
    } 
    else if (layer ==1 ) {
      prop_coord1        =      lt_prop_coord1_1[reducedAbsEta];
      prop_coord2        =      lt_prop_coord2_1[reducedAbsEta];
      res0_coord1        =      lt_res0_coord1_1[reducedAbsEta];
      res1_coord1        =      lt_res1_coord1_1[reducedAbsEta];
      res0_coord2        =      lt_res0_coord2_1[reducedAbsEta];
      res1_coord2        =      lt_res1_coord2_1[reducedAbsEta];
      res0_eta1          =      lt_res0_eta1_1[reducedAbsEta];
      res1_eta           =      lt_res1_eta_1[reducedAbsEta];
      res0_eta2          =      lt_res0_eta2_1[reducedAbsEta];
      valid              =      lt_prop_coord1_1_valid[reducedAbsEta];
      is_barrel          =      reducedAbsEta<barrelLimit1_ ? 1:0;

    }
    else if (layer ==2 ) {
      prop_coord1        =      lt_prop_coord1_2[reducedAbsEta];
      prop_coord2        =      lt_prop_coord2_2[reducedAbsEta];
      res0_coord1        =      lt_res0_coord1_2[reducedAbsEta];
      res1_coord1        =      lt_res1_coord1_2[reducedAbsEta];
      res0_coord2        =      lt_res0_coord2_2[reducedAbsEta];
      res1_coord2        =      lt_res1_coord2_2[reducedAbsEta];
      res0_eta1          =      lt_res0_eta1_2[reducedAbsEta];
      res1_eta           =      lt_res1_eta_2[reducedAbsEta];
      res0_eta2          =      lt_res0_eta2_2[reducedAbsEta];
      valid              =      lt_prop_coord1_2_valid[reducedAbsEta];
      is_barrel          =      reducedAbsEta<barrelLimit2_ ? 1:0;


    }
    else if (layer ==3 ) {
      prop_coord1        =      lt_prop_coord1_3[reducedAbsEta];
      prop_coord2        =      lt_prop_coord2_3[reducedAbsEta];
      res0_coord1        =      lt_res0_coord1_3[reducedAbsEta];
      res1_coord1        =      lt_res1_coord1_3[reducedAbsEta];
      res0_coord2        =      lt_res0_coord2_3[reducedAbsEta];
      res1_coord2        =      lt_res1_coord2_3[reducedAbsEta];
      res0_eta1          =      lt_res0_eta1_3[reducedAbsEta];
      res1_eta           =      lt_res1_eta_3[reducedAbsEta];
      res0_eta2          =      lt_res0_eta2_3[reducedAbsEta];
      valid              =      lt_prop_coord1_3_valid[reducedAbsEta];
      is_barrel          =      reducedAbsEta<barrelLimit3_ ? 1 :0;


    }
    else if (layer ==4 ) {
      prop_coord1        =      lt_prop_coord1_4[reducedAbsEta];
      prop_coord2        =      lt_prop_coord2_4[reducedAbsEta];
      res0_coord1        =      lt_res0_coord1_4[reducedAbsEta];
      res1_coord1        =      lt_res1_coord1_4[reducedAbsEta];
      res0_coord2        =      lt_res0_coord2_4[reducedAbsEta];
      res1_coord2        =      lt_res1_coord2_4[reducedAbsEta];
      res0_eta1          =      lt_res0_eta1_4[reducedAbsEta];
      res1_eta           =      lt_res1_eta_4[reducedAbsEta];
      res0_eta2          =      lt_res0_eta2_4[reducedAbsEta];
      valid              =      lt_prop_coord1_4_valid[reducedAbsEta];
      is_barrel          =      0;

    }

    propagation_t out;
    ap_int<BITSPHI> c1= track.phi()-prop_coord1*track.curvature()/1024;
    out.coord1 = c1/16;
    ap_int<BITSPHI> c2= prop_coord2*track.curvature()/1024;
    if (is_barrel==1)
      out.coord2 = -c2/16;
    else
      out.coord2 = (track.phi()-c2)/16;
    
    out.eta = track.eta()/(32);//to be updated with low pt data    
    out.valid=valid;
    
    ap_uint<BITSSTUBCOORD1-1> curvature2  = (track.curvature()*track.curvature())/(1<<19);
    out.sigma_coord1 = res0_coord1+res1_coord1*curvature2;
    if (out.sigma_coord1>128)
      out.sigma_coord1=128;
    out.sigma_coord2 = res0_coord2+res1_coord2*curvature2;
    if (out.sigma_coord2>128)
      out.sigma_coord2=128;

    out.sigma_eta1   = res0_eta1+res1_eta*curvature2;
    if (out.sigma_eta1>32)
      out.sigma_eta1=32;

    out.sigma_eta2   = res0_eta2+res1_eta*curvature2;
    if (out.sigma_eta2>32)
      out.sigma_eta2=32;
    out.is_barrel = is_barrel;

    if(verbose_==1)
      printf("Propagating to layer %d:is barrel=%d valid=%d coords=%d+-%d , %d +-%d etas = %d +- %d +-%d\n",int(layer),out.is_barrel.to_int(),valid.to_int(),out.coord1.to_int(),out.sigma_coord1.to_int(),out.coord2.to_int(),out.sigma_coord2.to_int(),out.eta.to_int(),out.sigma_eta1.to_int(),out.sigma_eta2.to_int());

    return out;
  }

  match_t match(const propagation_t prop,const l1t::MuonStubRef& stub) {

    if (verbose_==1) {
      printf("Matching to ");
      stub->print();
    }
  

    //Matching of Coord1 
    ap_uint<1> coord1Matched;
    ap_int<BITSSTUBCOORD1+1> deltaCoord1 = prop.coord1-stub->coord1();
    ap_uint<BITSSTUBCOORD1> absDeltaCoord1=0;
    if (deltaCoord1.sign())
      absDeltaCoord1=ap_uint<BITSSTUBCOORD1>(-deltaCoord1);
    else
      absDeltaCoord1=ap_uint<BITSSTUBCOORD1>(deltaCoord1);

    ap_int<BITSSTUBCOORD1+1> marginCoord1 = absDeltaCoord1-prop.sigma_coord1; 
    if (marginCoord1<=0 && (stub->quality() & 0x1) ) {
      coord1Matched=1;
    }
    else {
      coord1Matched=0;
    }

    if (verbose_==1)
      printf("Coord1 matched=%d delta=%d res=%d margin=%d\n",coord1Matched.to_int(),deltaCoord1.to_int(),prop.sigma_coord1.to_int(),marginCoord1.to_int());

    //Matching of Coord2 

    ap_uint<1> coord2Matched;
    ap_int<BITSSTUBCOORD2+1> deltaCoord2 = prop.coord2-stub->coord2();
    ap_uint<BITSSTUBCOORD2> absDeltaCoord2=0;
    if (deltaCoord2.sign())
      absDeltaCoord2=ap_uint<BITSSTUBCOORD2>(-deltaCoord2);
    else
      absDeltaCoord2=ap_uint<BITSSTUBCOORD2>(deltaCoord2);

    ap_int<BITSSTUBCOORD2+1> marginCoord2 = absDeltaCoord2-prop.sigma_coord2; 
    if (marginCoord2<=0 && (stub->quality() & 0x2) ) {
      coord2Matched=1;
    }
    else {
      coord2Matched=0;
    }

    if (verbose_==1)
      printf("Coord2 matched=%d delta=%d res=%d margin=%d\n",coord2Matched.to_int(),deltaCoord2.to_int(),prop.sigma_coord2.to_int(),marginCoord2.to_int());


    //Matching of Eta1



    ap_uint<1> eta1Matched;

    //if we have really bad quality[Barrel no eta]
    //increase the resolution
    ap_uint<BITSSTUBETA-1> prop_sigma_eta1 = stub->etaQuality()==0 ? ap_uint<BITSSTUBETA-1>(prop.sigma_eta1+200) : ap_uint<BITSSTUBETA-1>(prop_sigma_eta1);
    

    ap_int<BITSSTUBETA+1> deltaEta1 = prop.eta-stub->eta1();
    ap_uint<BITSSTUBETA> absDeltaEta1=0;
    if (deltaEta1.sign())
      absDeltaEta1=ap_uint<BITSSTUBETA>(-deltaEta1);
    else
      absDeltaEta1=ap_uint<BITSSTUBETA>(deltaEta1);

    ap_int<BITSSTUBETA+1> marginEta1 = absDeltaEta1-prop.sigma_eta1; 
    if (marginEta1<=0)
      eta1Matched=1;
    else
      eta1Matched=0;
    

    if (verbose_==1)
      printf("eta1 matched=%d delta=%d res=%d margin=%d\n",eta1Matched.to_int(),deltaEta1.to_int(),prop_sigma_eta1.to_int(),marginEta1.to_int());




    //Matching of Eta2

    ap_uint<1> eta2Matched;

    ap_int<BITSSTUBETA+1> deltaEta2 = prop.eta-stub->eta2();
    ap_uint<BITSSTUBETA> absDeltaEta2=0;
    if (deltaEta2.sign())
      absDeltaEta2=ap_uint<BITSSTUBETA>(-deltaEta2);
    else
      absDeltaEta2=ap_uint<BITSSTUBETA>(deltaEta2);

    ap_int<BITSSTUBETA+1> marginEta2 = absDeltaEta2-prop.sigma_eta2; 
    if (marginEta2<=0 && (stub->etaQuality() & 0x2))
      eta2Matched=1;
    else
      eta2Matched=0;
    
    match_t out;
    out.id=stub->id();


    if (verbose_==1)
      printf("eta2 matched=%d delta=%d res=%d margin=%d\n",eta2Matched.to_int(),deltaEta2.to_int(),prop.sigma_eta2.to_int(),marginEta2.to_int());

      
    //if barrel, coord1 has to always be matched, coord2 maybe and eta1 is needed if etaQ=0 or then the one that depends on eta quality 
    if (prop.is_barrel) {
      out.valid = (coord1Matched==1 && (eta1Matched==1||eta2Matched==1));
      if (out.valid==0) {
	out.quality=0;
      }
      else {
	out.quality = 64-deltaCoord1/2;
	if (coord2Matched==1)
	  out.quality=out.quality+64-deltaCoord2/2;
      }
    }
    //if endcap each coordinate is independent
    else {
      out.valid=(coord1Matched==1 &&eta1Matched==1) || (coord2Matched==1&&eta2Matched==1);
      if (out.valid==0) 
	out.quality=0;
      else {
	out.quality=0;
	if (coord1Matched==1)
	  out.quality+=64-deltaCoord1/2;
	if (coord2Matched==1)
	  out.quality+=64-deltaCoord2/2;
      } 

	
    }
    if (verbose_==1)
      printf("GlobalMatchQuality = %d\n",out.quality.to_int());
    out.stubRef = stub;
    return out;
  }


  match_t propagateAndMatch(const ConvertedTTTrack& track,const l1t::MuonStubRef& stub) {
    propagation_t prop = propagate(track,stub->tfLayer());
    return match(prop,stub);
  }
  

  match_t getBest(const std::vector<match_t> matches) {
    match_t best = matches[0];
    for (const auto& m : matches) {
      if (m.quality> best.quality)
	best=m;
    }

    return best;
  }
  

  PreTrackMatchedMuon processTrack(const ConvertedTTTrack& track , const std::vector<MuonROI> & rois) {
    std::vector<match_t> matchInfo0;
    std::vector<match_t> matchInfo1;
    std::vector<match_t> matchInfo2;
    std::vector<match_t> matchInfo3;
    std::vector<match_t> matchInfo4;


    if (verbose_==1 && rois.size()>0) {
      printf("-----------processing new track----------");
      track.print();
    }	     
    for (const auto& roi : rois)  {
      if (verbose_==1) {
	printf("New ROI with %d stubs \n",int(roi.stubs().size()));
      }
      for (const auto& stub: roi.stubs()) {
	match_t m = propagateAndMatch(track,stub);
	if (m.valid==1) {
	  if(stub->tfLayer()==0)
	    matchInfo0.push_back(m);
	  else if(stub->tfLayer()==1)
	    matchInfo1.push_back(m);
	  else if(stub->tfLayer()==2)
	    matchInfo2.push_back(m);
	  else if(stub->tfLayer()==3)
	    matchInfo3.push_back(m);
	  else if(stub->tfLayer()==4)
	    matchInfo4.push_back(m);
	  
	  if (roi.muonRef().isNonnull())
	    m.muRef = roi.muonRef();
	}
      }
    }
    
    uint quality=0;
    PreTrackMatchedMuon muon(track.charge(),track.pt(),track.eta(),track.phi(),track.z0(),track.d0());
    

    if (matchInfo0.size()>0) {
      match_t b = getBest(matchInfo0);
      if (b.valid) {
	muon.addStub(b.stubRef);
	if(!muon.muonRef().isNonnull())
	  muon.setMuonRef(b.muRef);
	quality+=b.quality;
      }
    }
    if (matchInfo1.size()>0) {
      match_t b = getBest(matchInfo1);
      if (b.valid) {
	muon.addStub(b.stubRef);
	if(!muon.muonRef().isNonnull())
	  muon.setMuonRef(b.muRef);
	quality+=b.quality;
      }
    }
    if (matchInfo2.size()>0) {
      match_t b = getBest(matchInfo2);
      if (b.valid) {
	muon.addStub(b.stubRef);
	if(!muon.muonRef().isNonnull())
	  muon.setMuonRef(b.muRef);
	quality+=b.quality;
      }
    }
    if (matchInfo3.size()>0) {
      match_t b = getBest(matchInfo3);
      if (b.valid) {
	muon.addStub(b.stubRef);
	if(!muon.muonRef().isNonnull())
	  muon.setMuonRef(b.muRef);
	quality+=b.quality;
      }
    }
    if (matchInfo4.size()>0) {

      match_t b = getBest(matchInfo4);
      if (b.valid) {
	muon.addStub(b.stubRef);
	if(!muon.muonRef().isNonnull())
	  muon.setMuonRef(b.muRef);
	quality+=b.quality;
      }
    }

      

    
    muon.setTrkPtr(track.trkPtr());
    muon.setQuality(quality);
    muon.setOfflineQuantities(track.offline_pt(),track.offline_eta(),track.offline_phi());
    return muon;
  }


  std::vector<PreTrackMatchedMuon> clean(std::vector<PreTrackMatchedMuon>&  muons) {
    std::vector<PreTrackMatchedMuon> out;
    if (muons.size()==0)
      return out;
    if (verbose_==1) {
      printf("-----Cleaning Up Muons in the same Nonant\n");
      printf("Before:\n");
    }
    for (uint i=0;i<muons.size();++i) {
      if (verbose_==1)
	muons[i].print();



      bool keep=true;
      for (uint j=0;j<muons.size();++j) {
	if (i==j)
	  continue;
	bool overlap=false;
	if (muons[i].stubID0()==muons[j].stubID0() && muons[i].stubID0()!=511)
	  overlap=true;
	if (muons[i].stubID1()==muons[j].stubID1() && muons[i].stubID1()!=511)
	  overlap=true;
	if (muons[i].stubID2()==muons[j].stubID2() && muons[i].stubID2()!=511)
	  overlap=true;
	if (muons[i].stubID3()==muons[j].stubID3() && muons[i].stubID3()!=511)
	  overlap=true;
	if (muons[i].stubID4()==muons[j].stubID4() && muons[i].stubID4()!=511)
	  overlap=true;

	if (overlap && muons[i].quality()<muons[j].quality())
	  keep=false;
      }
      if (keep) {
	if (verbose_==1)
	  printf("kept\n");
	out.push_back(muons[i]);
      }
      else {
	if (verbose_==1)
	  printf("discarded\n");
      }
    }
    return out;
  }





  };
}

#endif



