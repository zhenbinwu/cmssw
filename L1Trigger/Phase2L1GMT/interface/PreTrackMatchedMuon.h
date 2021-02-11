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
		      const int& d0
): 
    charge_(charge),
      pt_(pt),
      eta_(eta),
      phi_(phi),  
      z0_(z0),
      d0_(d0),
      quality_(0),

      stubID0_(511),
      stubID1_(511),
      stubID2_(511),
      stubID3_(511),
      stubID4_(511)
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
    void setQuality(uint quality) {
      quality_ = quality;
    }

    void setOfflineQuantities(float pt,float eta, float phi) {
      offline_pt_=pt;
      offline_eta_=eta;
      offline_phi_=phi;
    }

    void setMuonRef(const l1t::RegionalMuonCandRef& ref){
      muRef_=ref;
    }


    const l1t::RegionalMuonCandRef& muonRef() const {
      return muRef_;
    }
    void addStub(const l1t::MuonStubRef& stub){
      stubs_.push_back(stub);
      if (stub->tfLayer()==0)
	stubID0_ = stub->id();
      if (stub->tfLayer()==1)
	stubID1_ = stub->id();
      if (stub->tfLayer()==2)
	stubID2_ = stub->id();
      if (stub->tfLayer()==3)
	stubID3_ = stub->id();
      if (stub->tfLayer()==4)
	stubID4_ = stub->id();

    }

    const l1t::MuonStubRefVector& stubs() const {
      return stubs_;
    } 



    void setTrkPtr(const edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> >& trkPtr) {
      trkPtr_=trkPtr;
    }
    
    const edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> > trkPtr() const {
      return trkPtr_;
    }


    void print() const {
      printf("preconstructed muon  charge=%d pt=%f,%d eta=%f,%d phi=%f,%d z0=%d d0=%d quality=%d isGlobal=%d stubs: %d %d %d %d %d \n",charge_,offline_pt_,pt_,offline_eta_,eta_,offline_phi_,phi_,z0_,d0_,quality_,muRef_.isNonnull(),stubID0_,stubID1_,stubID2_,stubID3_,stubID4_);
    }




  private:
    uint          charge_;
    uint          pt_;
    int           eta_;
    int           phi_;
    int           z0_;
    int           d0_;
    uint          quality_;
    float         offline_pt_;
    float         offline_eta_;
    float         offline_phi_;
    int           stubID0_;
    int           stubID1_;
    int           stubID2_;
    int           stubID3_;
    int           stubID4_;
    l1t::MuonStubRefVector stubs_;
    l1t::RegionalMuonCandRef muRef_;
    edm::Ptr<TTTrack<Ref_Phase2TrackerDigi_> > trkPtr_;

  };
}

#endif
