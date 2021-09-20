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

#include "L1Trigger/Phase2L1GMT/interface/TopologicalAlgorithm.h"
#include "L1Trigger/Phase2L1GMT/interface/cos_dphi_LUT.h"
#include "L1Trigger/Phase2L1GMT/interface/cosh_deta_LUT.h"
#include "TLorentzVector.h"

namespace Phase2L1GMT {
  class Tauto3Mu : public TopoAlgo {
  public:
    Tauto3Mu(const edm::ParameterSet &iConfig);
    ~Tauto3Mu();
    Tauto3Mu(const Tauto3Mu &cpy);
    // Interface function
    bool GetTau3Mu(std::vector<l1t::TrackerMuon> &trkMus, std::vector<ConvertedTTTrack> &convertedTracks);

  private:
    bool Find3MuComb(std::vector<l1t::TrackerMuon> &trkMus);

    bool FindCloset3Mu(std::vector<std::pair<int, unsigned int> > &mu_phis,
                       std::vector<std::pair<unsigned, unsigned> > &nearby3mu);

    int Get3MuDphi(unsigned target, unsigned obj1, unsigned obj2);

    float Get3MuMass(unsigned target, unsigned obj1, unsigned obj2);

    float GetDiMass(const l1t::TrackerMuon &mu1, const l1t::TrackerMuon &mu2);
    float DumpInputs(unsigned target, unsigned obj1, unsigned obj2 );
    void DumpOutputs(float expMass, float emuMass);

    const int PT_MINV_LSB_RED = 2;
    const int pihalf = (1 << (BITSPHI -2));
    bool verbose_;
    bool dumpForHLS_;
    std::ofstream dumpOutput;
  };

  Tauto3Mu::Tauto3Mu(const edm::ParameterSet &iConfig):
    verbose_(iConfig.getParameter<int>("verbose")),
    dumpForHLS_(iConfig.getParameter<int>("IsodumpForHLS")) {
      dumpForHLS_ = true;
      if (dumpForHLS_) {
        dumpInput.open("ta3mu_input_data.txt", std::ofstream::out);
        dumpOutput.open("HLS_m3mu_testebench.txt", std::ofstream::out);
      }
    }

  Tauto3Mu::~Tauto3Mu() {
    if (dumpForHLS_) {
      dumpInput.close();
      dumpOutput.close();
    }
  }

  Tauto3Mu::Tauto3Mu(const Tauto3Mu &cpy) {}


  float Tauto3Mu::DumpInputs(unsigned target, unsigned obj1, unsigned obj2 ) {
    const float muMass = 0.105;

    l1t::TrackerMuon &mu0 = trkMus->at(target);
    l1t::TrackerMuon &mu1 = trkMus->at(obj1);
    l1t::TrackerMuon &mu2 = trkMus->at(obj2);
    TLorentzVector v0(0, 0, 0, 0);
    TLorentzVector v1(0, 0, 0, 0);
    TLorentzVector v2(0, 0, 0, 0);
    v0.SetPtEtaPhiM(mu0.hwPt() * LSBpt, mu0.hwEta()*LSBeta, mu0.hwPhi()*LSBphi, 0);
    v1.SetPtEtaPhiM(mu1.hwPt() * LSBpt, mu1.hwEta()*LSBeta, mu1.hwPhi()*LSBphi, 0);
    v2.SetPtEtaPhiM(mu2.hwPt() * LSBpt, mu2.hwEta()*LSBeta, mu2.hwPhi()*LSBphi, 0);
    TLorentzVector minv_massless = v0+v1+v2;

    v0.SetPtEtaPhiM(mu0.hwPt() * LSBpt, mu0.hwEta()*LSBeta, mu0.hwPhi()*LSBphi, muMass);
    v1.SetPtEtaPhiM(mu1.hwPt() * LSBpt, mu1.hwEta()*LSBeta, mu1.hwPhi()*LSBphi, muMass);
    v2.SetPtEtaPhiM(mu2.hwPt() * LSBpt, mu2.hwEta()*LSBeta, mu2.hwPhi()*LSBphi, muMass);
    TLorentzVector minv = v0+v1+v2;

    dumpInput
      << v0.Pt() <<" 0 " << v0.Phi()<<" "<< v0.Eta() <<" "
      << v1.Pt() <<" 0 " << v1.Phi()<<" "<< v1.Eta() <<" "
      << v2.Pt() <<" 0 " << v2.Phi()<<" "<< v2.Eta() <<" "
      << minv.M() <<" "<< minv_massless.M()<<std::endl;
    return minv_massless.M();

  }

  void Tauto3Mu::DumpOutputs(float expMass, float emuMass) {
        dumpOutput <<  emuMass<<" 0 "<< expMass <<" 0 "<< endl;
  }

  // ===  FUNCTION  ============================================================
  //         Name:  Tauto3Mu::GetTau3Mu
  //  Description:
  // ===========================================================================
  bool Tauto3Mu::GetTau3Mu(std::vector<l1t::TrackerMuon> &trkMus, std::vector<ConvertedTTTrack> &convertedTracks) {
    if (trkMus.size() < 3)
      return false;

    load(trkMus, convertedTracks);
    Find3MuComb(trkMus);
    return true;
  }  // -----  end of function Tauto3Mu::GetTau3Mu  -----

  // ===  FUNCTION  ============================================================
  //         Name:  Tauto3Mu::Find3MuComb
  //  Description:
  // ===========================================================================
  bool Tauto3Mu::Find3MuComb(std::vector<l1t::TrackerMuon> &trkMus) {
    // vector< phi, index of trackerMuon >
    std::vector<std::pair<int, unsigned int> > mu_phis;
    for (unsigned i = 0; i < trkMus.size(); ++i) {
      mu_phis.push_back(std::make_pair(trkMus.at(i).hwPhi(), i));
    }

    std::sort(mu_phis.begin(), mu_phis.end());

    std::vector<std::pair<unsigned, unsigned> > nearby3mu;
    std::vector<int> mu3mass;
    FindCloset3Mu(mu_phis, nearby3mu);

    for (unsigned i = 0; i < trkMus.size(); ++i) {
      float trimass = Get3MuMass(i, nearby3mu.at(i).first, nearby3mu.at(i).second);
      mu3mass.push_back(trimass);

      if (dumpForHLS_) {
        float expMass = DumpInputs(i, nearby3mu.at(i).first, nearby3mu.at(i).second);
        float cvtmass = sqrt(trimass) * LSBpt *(1<<PT_MINV_LSB_RED);
        //float cvtmass = trimass;
        DumpOutputs(expMass, cvtmass);
      }
    }

    return true;
  }  // -----  end of function Tauto3Mu::Find3MuComb  -----

  // ===  FUNCTION  ============================================================
  //         Name:  Tauto3Mu::Get3MuMass
  //  Description:
  // ===========================================================================
  float Tauto3Mu::Get3MuMass(unsigned target, unsigned obj1, unsigned obj2) {
    float mass12 = GetDiMass(trkMus->at(target), trkMus->at(obj1));
    float mass23 = GetDiMass(trkMus->at(obj1), trkMus->at(obj2));
    float mass31 = GetDiMass(trkMus->at(obj2), trkMus->at(target));

    return mass12 + mass23 + mass31;
    //typedef ap_ufixed<45, 26, AP_TRN, AP_SAT> hw_minv2over2_t; // m^2/2, LSB : 25 MeV^2 x (2^2)^2
    //hw_minv2over2_t temp(mass12);
    //return temp;
  }  // -----  end of function Tauto3Mu::Get3MuMass  -----

  // ===  FUNCTION  ============================================================
  //         Name:  Tauto3Mu::GetDiMass
  //  Description:
  // ===========================================================================
  float Tauto3Mu::GetDiMass(const l1t::TrackerMuon &mu1, const l1t::TrackerMuon &mu2) {
    int ptprod = (mu1.hwPt() >> PT_MINV_LSB_RED) * (mu2.hwPt()  >> PT_MINV_LSB_RED);

    int deta = deltaEta(mu1.hwEta(), mu2.hwEta());
    float coshdeta = cosh_deta_lut[deta>>3];

    int dphi = deltaPhi(mu1.hwPhi(), mu2.hwPhi());
    bool abovepihalf = dphi > pihalf;
    if (abovepihalf)
        dphi = dphi - pihalf;
    float cosdphi = cos_dphi_lut[dphi];
    if (abovepihalf)
        cosdphi = -1 * cosdphi;

    float mass = 2 * ptprod * ( coshdeta - cosdphi);
    return mass;
  }  // -----  end of function Tauto3Mu::GetDiMass  -----

  // ===  FUNCTION  ============================================================
  //         Name:  Tauto3Mu::FindCloset3Mu
  //  Description:
  // ===========================================================================
  bool Tauto3Mu::FindCloset3Mu(std::vector<std::pair<int, unsigned int> > &mu_phis,
                               std::vector<std::pair<unsigned, unsigned> > &nearby3mu) {
    nearby3mu.clear();

    std::vector<std::pair<int, unsigned int> > temp(mu_phis);

    // Round the last 2 to first element of vector
    temp.insert(temp.begin(), mu_phis.back());
    temp.insert(temp.begin(), *(mu_phis.rbegin() + 1));
    // Append the first two element to vector
    temp.push_back(mu_phis.front());
    temp.push_back(*(mu_phis.begin() + 1));

    for (unsigned i = 2; i < temp.size() - 2; ++i) {
      int combleft = Get3MuDphi(temp[i].second, temp[i - 1].second, temp[i - 2].second);
      std::pair<unsigned, unsigned> neighbors(temp[i - 1].second, temp[i - 2].second);
      int mincomb(combleft);

      int combcenter = Get3MuDphi(temp[i].second, temp[i - 1].second, temp[i + 1].second);
      if (combcenter < mincomb) {
        neighbors = std::make_pair(temp[i - 1].second, temp[i + 1].second);
        mincomb = combcenter;
      }

      int combright = Get3MuDphi(temp[i].second, temp[i + 1].second, temp[i + 2].second);
      if (combright < mincomb) {
        neighbors = std::make_pair(temp[i + 1].second, temp[i + 2].second);
      }

      nearby3mu.push_back(neighbors);
    }

    return true;
  }  // -----  end of function Tauto3Mu::FindCloset3Mu  -----

  int Tauto3Mu::Get3MuDphi(unsigned target, unsigned obj1, unsigned obj2) {
    int dPhi1 = deltaPhi(trkMus->at(target).hwPhi(), trkMus->at(obj1).hwPhi());
    int dPhi2 = deltaPhi(trkMus->at(target).hwPhi(), trkMus->at(obj2).hwPhi());
    return dPhi1 + dPhi2;
  }
}  // namespace Phase2L1GMT

#endif  // ----- #ifndef PHASE2GMT_TAUTO3MU -----
