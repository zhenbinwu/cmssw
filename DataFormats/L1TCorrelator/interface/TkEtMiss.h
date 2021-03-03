#ifndef DataFormatsL1TCorrelator_TkEtMiss_h
#define DataFormatsL1TCorrelator_TkEtMiss_h

#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/L1TCorrelator/interface/TkPrimaryVertex.h"

namespace l1t {
  class TkEtMiss : public L1Candidate {
  public:
    enum EtMissType { kMET, kMHT, kNumTypes,hwMET };
    TkEtMiss();
    TkEtMiss(const LorentzVector& p4,
             EtMissType type,
             const double& etTotal,
             const double& etMissPU,
             const double& etTotalPU,
             const edm::Ref<TkPrimaryVertexCollection>& aVtxRef = edm::Ref<TkPrimaryVertexCollection>(),
             int bx = 0);

    TkEtMiss(const LorentzVector& p4,
             EtMissType type,
             const double& etTotal,
             const double& etMissPU,
             const double& etTotalPU,
             int bx = 0);

    TkEtMiss(const LorentzVector& p4,
             EtMissType type,
             const unsigned int& Etmiss,
             const unsigned int& EtPhi,
             const unsigned int& NumTracks ,
             const double& hwEtScale,
             const double& hwPhiScale,
             int bx = 0);

    TkEtMiss(const LorentzVector& p4,
             EtMissType type,
             const double& Etmiss,
             const double& EtPhi,
             const unsigned int& NumTracks ,
             int bx = 0);

    // ---------- const member functions ---------------------
    EtMissType type() const { return type_; }  // kMET or kMHT
    // For type = kMET, this is |MET|; for type = kMHT, this is |MHT|
    double etMiss() const { return et(); }
    // For type = kMET, this is total ET; for type = kMHT, this is total HT
    double etTotal() const { return etTot_; }
    // EtMiss and EtTot from PU vertices
    double etMissPU() const { return etMissPU_; }
    double etTotalPU() const { return etTotalPU_; }
    int bx() const { return bx_; }
    const edm::Ref<TkPrimaryVertexCollection>& vtxRef() const { return vtxRef_; }


    double hwEtmiss() const {return hwEtmiss_;}
    double hwEtPhi()  const {return hwEtphi_;}
    unsigned int hwNumTracks() const {return hwNumTracks_;}

    unsigned int hwEtmissBits() const {return hwiEtmiss_;}
    unsigned int hwEtPhiBits()  const {return hwiEtphi_;}

    // ---------- member functions ---------------------------
    void setEtTotal(const double& etTotal) { etTot_ = etTotal; }
    void setBx(int bx) { bx_ = bx; }

    void ihwtofloat();
    void hwfloattoi();


  private:
    // ---------- member data --------------------------------
    EtMissType type_;
    double etTot_;
    double etMissPU_;
    double etTotalPU_;
    edm::Ref<TkPrimaryVertexCollection> vtxRef_;
    


    unsigned int hwiEtmiss_;
    unsigned int hwiEtphi_;
    double hwEtmiss_;
    double hwEtphi_;
    unsigned int hwNumTracks_;
    

    double hwEtScale_;
    double hwPhiScale_;

    int bx_;

  };
}  // namespace l1t

#endif
