#ifndef PHASE2GMT_TEMPORARY_ASSOCIATOR
#define PHASE2GMT_TEMPORARY_ASSOCIATOR

#include "ap_int.h"
#include "L1Trigger/Phase2L1GMT/interface/MuonROI.h"
#include "DataFormats/L1TMuonPhase2/interface/L1TPhase2GMTStub.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

namespace Phase2L1GMT {

  class ROITempAssociator {
  public:
  ROITempAssociator(const edm::ParameterSet& iConfig)  {}
    ~ROITempAssociator() {}

    std::vector<MuonROI> associate(int bx,const l1t::ObjectRefBxCollection<l1t::RegionalMuonCand>& muons,const  L1TPhase2GMTStubRefVector& stubs )
      {
	std::vector<MuonROI> out;
	L1TPhase2GMTStubRefVector usedStubs;


	
	if (muons.size()>0) {
	  for (unsigned int i=0;i<muons.size(bx);++i) {
	    const l1t::RegionalMuonCandRef& mu = muons.at(bx,i);
	    uint pt = mu->hwPt();
	    uint charge = mu->hwSign();

	    float eta = mu->hwEta()*0.010875;

	    int globalPhi=0;
	    if (mu->trackFinderType()==l1t::bmtf) {
	      globalPhi = mu->processor()*48+mu->hwPhi()-24;
	    }
	    else {
	      globalPhi = mu->processor()*96+mu->hwPhi()+24;
	    }

	    float phi = globalPhi*2*M_PI/576.0;
	    if (phi>(M_PI))
	      phi=phi-2*M_PI;
	    if (phi<(-M_PI))
	      phi=phi+2*M_PI;


	    MuonROI roi(bx,charge,pt,1);
	    L1TPhase2GMTStubRefVector cleanedStubs = clean(stubs,usedStubs);
	  
	    for (unsigned int layer=0;layer<=4;++layer) {
	      L1TPhase2GMTStubRefVector selectedStubs = selectLayerBX(cleanedStubs,bx,layer);
	      int bestStubINT=-1;
	      float dPhi=1000.0;

	      for (uint i=0;i<selectedStubs.size();++i) {
		const L1TPhase2GMTStubRef& stub = selectedStubs[i];
		float deltaPhi = (stub->quality() & 0x1) ? stub->offline_coord1()-phi :stub->offline_coord2()-phi ;
		if (deltaPhi>M_PI)
		  deltaPhi=deltaPhi-2*M_PI;
		if (deltaPhi<-M_PI)
		  deltaPhi=deltaPhi+2*M_PI;
		deltaPhi=fabs(deltaPhi);
		float deltaEta= (stub->etaQuality() & 0x1) ?  fabs(stub->offline_eta1()-eta) : fabs(stub->offline_eta2()-eta);
		if (deltaPhi<(M_PI/6.0) && deltaEta<0.3 && deltaPhi<dPhi) {
		  dPhi=deltaPhi;
		  bestStubINT=i;
		}
	      }
	      if (bestStubINT>=0) {
		roi.addStub(selectedStubs[bestStubINT]);
		usedStubs.push_back(selectedStubs[bestStubINT]);
	      }
	    }
	    out.push_back(roi);
	  }
	}
	//Now the stubs only . Find per layer

	L1TPhase2GMTStubRefVector cleanedStubs = clean(stubs,usedStubs);
	

	
	while (cleanedStubs.size()>0) {
	  MuonROI roi(bx,0,0,0);
	  roi.addStub(cleanedStubs[0]);
	  usedStubs.push_back(cleanedStubs[0]);
	  for (unsigned int layer=0;layer<=4;++layer) {
	    if (layer==cleanedStubs[0]->tfLayer())
	      continue;
	    L1TPhase2GMTStubRefVector selectedStubs = selectLayerBX(cleanedStubs,bx,layer);
	    if (selectedStubs.size()!=0) {
	      roi.addStub(selectedStubs[0]);
	      usedStubs.push_back(selectedStubs[0]);
	    }
	  }
	  if (roi.stubs().size()>0) 
	    out.push_back(roi);
	  cleanedStubs = clean(cleanedStubs,usedStubs);
	}



	return out;
      }

  private:

    L1TPhase2GMTStubRefVector selectLayerBX(const L1TPhase2GMTStubRefVector& all, int bx, uint layer) {
      L1TPhase2GMTStubRefVector out;
      for (const auto& stub: all) {
	if (stub->bxNum()==bx && stub->tfLayer()==layer)
	  out.push_back(stub);
      }
      return out;
    }

    L1TPhase2GMTStubRefVector clean(const L1TPhase2GMTStubRefVector& all, const L1TPhase2GMTStubRefVector& used) {
      L1TPhase2GMTStubRefVector out;
      for (const auto& stub : all) {
	bool keep=true;
	for (const auto& st : used)  {
	  if (st==stub)  {
	    keep=false;
	    break;
	  }
	}
	if (keep)
	  out.push_back(stub);
      }
      return out;
    }


    



  };
}

#endif



