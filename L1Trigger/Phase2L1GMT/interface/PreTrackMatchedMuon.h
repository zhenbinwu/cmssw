#ifndef PHASE2GMT_PRETRACKMATCHEDMUON
#define PHASE2GMT_PRETRACKMATCHEDMUON



namespace Phase2L1GMT {

  class PreTrackMatchedMuon {

  public:
  PreTrackMatchedMuon(const uint& charge,
		      const uint & pt,
		      const int& eta,
		      const int& phi,
		      const int& z0,
		      const int& d0,
		      const int& quality,
		      const uint& isGlobal,
		      const int& stubID0,   
		      const int& stubID1,   
		      const int& stubID2,   
		      const int& stubID3,   
		      const int& stubID4   
): 
    charge_(charge),
      pt_(pt),
      eta_(eta),
      phi_(phi),  
      z0_(z0),
      d0_(d0),
      quality_(quality),
      isGlobal_(isGlobal),
      stubID0_(stubID0),
      stubID1_(stubID1),
      stubID2_(stubID2),
      stubID3_(stubID3),
      stubID4_(stubID4)
    {
    }

    const uint charge() const {
      return charge_;
    }
    const uint pt() const {
      return pt_;
    }
    const int eta() const {
      return eta_;
    }
    const int phi() const {
      return phi_;
    }
    const int z0() const {
      return z0_;
    }
    const int d0() const {
      return d0_;
    }
    const int quality() const {
      return quality_;
    }
    const int offline_pt() const {
      return offline_pt_;
    }
    const float offline_eta() const {
      return offline_eta_;
    }
    const float offline_phi() const {
      return offline_phi_;
    }

    const int stubID0() const {
      return stubID0_;
    }
    const int stubID1() const {
      return stubID1_;
    }
    const int stubID2() const {
      return stubID2_;
    }
    const int stubID3() const {
      return stubID3_;
    }
    const int stubID4() const {
      return stubID4_;
    }

    void setOfflineQuantities(float pt,float eta, float phi) {
      offline_pt_=pt;
      offline_eta_=eta;
      offline_phi_=phi;
    }

    void setMuonReference(const l1t::RegionalMuonCandRef& ref){
      muRef_=ref;
    }


    const l1t::RegionalMuonCandRef& muonRef() const {
      return muRef_;
    }
    void addStub(const L1TPhase2GMTStubRef& stub){
      stubs_.push_back(stub);
    }

    const L1TPhase2GMTStubRefVector& stubs() const {
      return stubs_;
    } 

    void print() const {
      printf("converted track charge=%d pt=%f,%d eta=%f,%d phi=%f,%d z0=%d d0=%d quality=%d\n",charge_,offline_pt_,pt_,offline_eta_,eta_,offline_phi_,phi_,z0_,d0_,quality_);
    }




  private:
    uint          charge_;
    uint          pt_;
    int           eta_;
    int           phi_;
    int           z0_;
    int           d0_;
    uint          quality_;
    uint          isGlobal_;
    float         offline_pt_;
    float         offline_eta_;
    float         offline_phi_;
    int           stubID0_;
    int           stubID1_;
    int           stubID2_;
    int           stubID3_;
    int           stubID4_;
    L1TPhase2GMTStubRefVector stubs_;
    l1t::RegionalMuonCandRef muRef_;

  };
}

#endif
