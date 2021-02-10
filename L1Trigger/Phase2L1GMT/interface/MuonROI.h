#ifndef PHASE2GMT_MUONROI
#define PHASE2GMT_MUONROI
#include <iosfwd>
#include "L1Trigger/Phase2L1GMT/interface/Constants.h"
#include "DataFormats/L1TMuonPhase2/interface/L1TPhase2GMTStub.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1Trigger/interface/L1TObjComparison.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

namespace Phase2L1GMT {

  class MuonROI {

  public:
  MuonROI(int bx,
	  uint charge,
	  uint pt,
	  uint quality
	  ): 
    bx_(bx),
    charge_(charge),
    pt_(pt),
    quality_(quality)
    {
    }

    const int bx() const {
      return bx_;
    }


    const uint charge() const {
      return charge_;
    }

    const uint pt() const {
      return pt_;
    }
    const int quality() const {
      return quality_;
    }

    const float offline_pt() const {
      return offline_pt_;
    }

    void setOfflinePt(float pt) {
      offline_pt_=pt;
    }

    void addStub(const L1TPhase2GMTStubRef& stub){
      stubs_.push_back(stub);
    }

    void setMuonReference(const l1t::RegionalMuonCandRef& ref){
      muRef_=ref;
    }


    const l1t::RegionalMuonCandRef& muonRef() const {
      return muRef_;
    }

    friend std::ostream& operator<<(std::ostream& s, const MuonROI& id) {
      s.setf(ios::right,ios::adjustfield);
      s << "ROI:" << " "
        << "BX: "              << setw(5) << id.bx_        << " "
        << "charge:"           << setw(5) << id.charge_    << " "
        << "pt:"               << setw(5) << id.pt_        << " "
        << "quality:"          << setw(5) << id.quality_    << " "
        << "offline pt:"       << setw(5) << id.offline_pt_;
      return s;
    }

    const L1TPhase2GMTStubRefVector& stubs() const {
      return stubs_;
    } 
  private:
    int                       bx_;
    uint                      charge_;
    uint                      pt_;
    uint                      quality_;
    float                     offline_pt_;
    L1TPhase2GMTStubRefVector stubs_;
    l1t::RegionalMuonCandRef muRef_;

  };
}

#endif
