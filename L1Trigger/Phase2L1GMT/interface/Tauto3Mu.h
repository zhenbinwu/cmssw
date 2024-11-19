// ===========================================================================
//
//       Filename:  Tauto3Mu.h
//
//    Description:
//
//        Version:  1.0
//        Created:  03/15/2021 07:33:59 PM
//       Revision:  none
//       Compiler:  g++
//
//         Author:  Zhenbin Wu, zhenbin.wu@gmail.com
//
// ===========================================================================

#ifndef PHASE2GMT_TAUTO3MU
#define PHASE2GMT_TAUTO3MU

#include "cos_dphi_LUT.h"
#include "cosh_deta_LUT.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/L1TMuonPhase2/interface/Tau23Mu.h"
#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "DataFormats/L1TMuonPhase2/interface/Constants.h"
#include <fstream>

namespace Phase2L1GMT {
  //class Tauto3Mu : public TopoAlgo {
  class Tauto3Mu {
  public:
    Tauto3Mu(const edm::ParameterSet &iConfig);
    ~Tauto3Mu();
    Tauto3Mu(const Tauto3Mu &cpy);
    // Interface function
    std::vector<l1t::Tau23Mu> GetTau3Mu(const std::vector<l1t::TrackerMuon> &trkMus);
    inline int deltaEta(const int eta1, const int eta2);
    inline int deltaPhi(int phi1, int phi2);

  private:
    bool Find3MuComb(std::vector<l1t::TrackerMuon> &trkMus);

    bool FindCloset3Mu(std::vector<std::vector<unsigned> > &nearby3mu);

    int Get3MuDphi(unsigned target, unsigned obj1, unsigned obj2);

    float Get3MuMass(unsigned target, unsigned obj1, unsigned obj2);

    float GetDiMass(const l1t::TrackerMuon &mu1, const l1t::TrackerMuon &mu2);
    float DumpInputs(unsigned target, unsigned obj1, unsigned obj2);

    void DumpOutputs(int Nevt, unsigned obj0, unsigned obj1, unsigned obj2, float expMass, float emuMass);
    void DumpOutputs(unsigned obj0, unsigned obj1, unsigned obj2, float expMass, float emuMass);

    const int PT_MINV_LSB_RED = 2;
    const int pihalf = (1 << (BITSPHI - 2));
    int Nevt = -1;
    bool verbose_;
    int hwdRcut_;
    int hwdzcut_;
    bool dodRcut_;
    bool dodzcut_;
    bool dumpForHLS_;
    std::ofstream dumpOutput;
    std::vector<l1t::Tau23Mu> outtau;
    std::ofstream dumpInput;
    std::vector<l1t::TrackerMuon> *trkMus;
  };

}  // namespace Phase2L1GMT

#endif  // ----- #ifndef PHASE2GMT_TAUTO3MU -----
