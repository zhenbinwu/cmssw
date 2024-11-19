// ===========================================================================
// 
//       Filename:  Tau23Mu.h
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  09/22/2021 03:46:35 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Zhenbin Wu (benwu), zhenbin.wu@gmail.com
// 
// ===========================================================================


#ifndef dataformatsl1tmuonphase2_tau23mu_h
#define dataformatsl1tmuonphase2_tau23mu_h

#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

namespace l1t {

  class Tau23Mu;

  typedef std::vector<Tau23Mu> Tau23MuCollection;
  typedef edm::Ref<Tau23MuCollection> Tau23MuRef;
  typedef std::vector<edm::Ref<Tau23MuCollection> > Tau23MuRefVector;

  class Tau23Mu : public L1Candidate {
  public:
    Tau23Mu();

    Tau23Mu(uint pt, int eta, int phi, uint mass=0, uint quality=0, int mu1=0, int mu2=0, int mu3=0);

    ~Tau23Mu() override;

    void setWord(uint64_t word) { word_ = word; }
    void setHwMu1(int mu1) { mu1_ = mu1; }
    void setHwMu2(int mu2) { mu2_ = mu2; }
    void setHwMu3(int mu3) { mu3_ = mu3; }
    void setHwMass(int mass) { hwMass_ = mass; }

    const uint64_t word() const { return word_; }
    const uint64_t hwMass() const { return hwMass_; }
    const int hwMu1() const { return mu1_; }
    const int hwMu2() const { return mu2_; }
    const int hwMu3() const { return mu3_; }

    bool operator<(const Tau23Mu& other) const {
      if (hwPt() == other.hwPt())
        return (hwEta() < other.hwEta());
      else
        return (hwPt() < other.hwPt());
    }
    bool operator>(const Tau23Mu& other) const {
      if (hwPt() == other.hwPt())
        return (hwEta() > other.hwEta());
      else
        return (hwPt() > other.hwPt());
    }

    void print() const;

  private:
    int hwMass_;
    int mu1_;
    int mu2_;
    int mu3_;
    uint64_t word_;
  };
}  // namespace l1t

#endif   // ----- #ifndef dataformatsl1tmuonphase2_tau23mu_h -----
