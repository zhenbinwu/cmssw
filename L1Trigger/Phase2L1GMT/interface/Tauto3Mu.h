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
    std::vector<l1t::Tau23Mu> GetTau3Mu(std::vector<l1t::TrackerMuon> &trkMus, std::vector<ConvertedTTTrack> &convertedTracks);

  private:
    bool Find3MuComb(std::vector<l1t::TrackerMuon> &trkMus);

    bool FindCloset3Mu(std::vector<std::vector<unsigned> > &nearby3mu);

    int Get3MuDphi(unsigned target, unsigned obj1, unsigned obj2);

    float Get3MuMass(unsigned target, unsigned obj1, unsigned obj2);

    float GetDiMass(const l1t::TrackerMuon &mu1, const l1t::TrackerMuon &mu2);
    float DumpInputs(unsigned target, unsigned obj1, unsigned obj2 );

    void DumpOutputs(int Nevt, unsigned obj0, unsigned obj1, unsigned obj2, float expMass, float emuMass);
    void DumpOutputs(unsigned obj0, unsigned obj1, unsigned obj2, float expMass, float emuMass);

    const int PT_MINV_LSB_RED = 2;
    const int pihalf = (1 << (BITSPHI -2));
    int Nevt=-1;
    bool verbose_;
    bool dumpForHLS_;
    std::ofstream dumpOutput;
    std::vector<l1t::Tau23Mu> outtau;
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

  void Tauto3Mu::DumpOutputs(unsigned obj0, unsigned obj1, unsigned obj2, float expMass, float emuMass) {
    dumpOutput <<  
      trkMus->at(obj0).hwPt() << " "<< trkMus->at(obj0).hwPhi() << " "<< trkMus->at(obj0).hwEta() << " "<<
      trkMus->at(obj1).hwPt() << " "<< trkMus->at(obj1).hwPhi() << " "<< trkMus->at(obj1).hwEta() << " "<<
      trkMus->at(obj2).hwPt() << " "<< trkMus->at(obj2).hwPhi() << " "<< trkMus->at(obj2).hwEta() << " "<<
      emuMass<<"  "<< expMass << endl;
  }

  void Tauto3Mu::DumpOutputs(int Nevt, unsigned obj0, unsigned obj1, unsigned obj2, float expMass, float emuMass) {
    dumpOutput << Nevt <<" : ";
    DumpOutputs(obj0, obj1, obj2, expMass, emuMass);
  }

  // ===  FUNCTION  ============================================================
  //         Name:  Tauto3Mu::GetTau3Mu
  //  Description:
  // ===========================================================================
  std::vector<l1t::Tau23Mu> Tauto3Mu::GetTau3Mu(std::vector<l1t::TrackerMuon> &trkMus, std::vector<ConvertedTTTrack> &convertedTracks) {
    outtau.clear();
    Nevt++;

    if (trkMus.size() < 3)
      return outtau;

    load(trkMus, convertedTracks);
    Find3MuComb(trkMus);
    return outtau;
  }  // -----  end of function Tauto3Mu::GetTau3Mu  -----

  // ===  FUNCTION  ============================================================
  //         Name:  Tauto3Mu::Find3MuComb
  //  Description:
  // ===========================================================================
  bool Tauto3Mu::Find3MuComb(std::vector<l1t::TrackerMuon> &trkMus) {
    std::vector<std::vector<unsigned> > nearby3mu;
    FindCloset3Mu(nearby3mu);

    for (unsigned i = 0; i < trkMus.size(); ++i) {

      l1t::TrackerMuon &seedmu = trkMus[nearby3mu[i][0]];
      float trimass = Get3MuMass(nearby3mu[i][0], nearby3mu[i][1], nearby3mu[i][2]);
      int outmass(int(trimass) >> 3);

      l1t::Tau23Mu temp(seedmu.hwPt(), seedmu.hwEta(), seedmu.hwPhi(), trimass, 
          0, nearby3mu[i][0], nearby3mu[i][1], nearby3mu[i][2]);

      outtau.push_back(temp);

      if (dumpForHLS_) {
        float expMass = DumpInputs(nearby3mu[i][0], nearby3mu[i][1], nearby3mu[i][2]);
        float cvtmass = sqrt(2*trimass) * LSBpt *(1<<PT_MINV_LSB_RED);
        //std::cout << " expMass " << expMass <<" cvtMass " << cvtmass << " outmass " << outmass <<" expIntMAss " << cvtmass*cvtmass/0.25 <<" diff " << cvtmass*cvtmass/outmass << std::endl;
        DumpOutputs(Nevt, nearby3mu[i][0], nearby3mu[i][1], nearby3mu[i][2], expMass, cvtmass);
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

    ap_ufixed<45, 26, AP_TRN, AP_SAT> temp(mass12 + mass23 + mass31);
    if (verbose_)
    {
      std::cout <<  target <<" "<<obj1 <<" "<<obj2 << std::endl;
      std::cout << "mass12 " <<  mass12 <<", "
        << "mass23 " <<  mass23 <<", "
        << "mass31 " <<  mass31 <<", "
        <<" total " << temp <<", "
        <<" float " << sqrt(2*temp.to_float()) * LSBpt *(1<<PT_MINV_LSB_RED)
        << std::endl;
    }
    return temp;
  }  // -----  end of function Tauto3Mu::Get3MuMass  -----

  // ===  FUNCTION  ============================================================
  //         Name:  Tauto3Mu::GetDiMass
  //  Description:
  // ===========================================================================
  float Tauto3Mu::GetDiMass(const l1t::TrackerMuon &mu1, const l1t::TrackerMuon &mu2) {
    int ptprod = (mu1.hwPt() >> PT_MINV_LSB_RED) * (mu2.hwPt()  >> PT_MINV_LSB_RED);

    int deta = deltaEta(mu1.hwEta(), mu2.hwEta());
    float coshdeta = cosh_deta_lut[deta>>3].to_float();

    int dphi = deltaPhi(mu1.hwPhi(), mu2.hwPhi());
    bool abovepihalf = dphi > pihalf;
    if (abovepihalf)
        dphi = dphi - pihalf;
    float cosdphi = cos_dphi_lut[dphi].to_float();
    if (abovepihalf)
        cosdphi = -1 * cosdphi;

    float mass = ptprod * ( coshdeta - cosdphi);
    if (verbose_)
    {
         std::cout << " ..... "
                   << " pt1 "        << mu1.hwPt()
                   << " pt2 "        << mu2.hwPt()
                   << " ptprod "     << ptprod
                   << "; deta "      << deta
                   << " coshdeta "   << coshdeta
                   << "; phi1 "      << mu1.hwPhi()
                   << " phi2 "       << mu2.hwPhi()
                   << " dphi "       << dphi
                   << " cosdphi "    << cosdphi
                   << " ; ang_diff " << ( coshdeta - cosdphi)
                   << " minv2o2 "    << mass
                   << std::endl;
    }
    return mass;
  }  // -----  end of function Tauto3Mu::GetDiMass  -----

  // ===  FUNCTION  ============================================================
  //         Name:  Tauto3Mu::FindCloset3Mu
  //  Description:
  // ===========================================================================
  bool Tauto3Mu::FindCloset3Mu(std::vector<std::vector<unsigned> > &nearby3mu) {

    nearby3mu.clear();


    for (unsigned i = 0; i < trkMus->size(); ++i)
    {
      std::pair<unsigned, unsigned > minComp;
      int mindR = 99999;

      for (unsigned j = 0; j < trkMus->size(); ++j)
      {
        if (i == j) continue;
        for (unsigned k = 0; k < trkMus->size(); ++k)
        {
          if (i == k) continue;
          if (j == k) continue;
          int comp = Get3MuDphi(i, j, k);
          if (comp < mindR)
          {
            mindR = comp;
            minComp=std::make_pair(j, k);
          }
        }
      }
      std::vector<unsigned> tempv={i, minComp.first, minComp.second};
      nearby3mu.push_back(tempv);
    }
    return true;
  }  // -----  end of function Tauto3Mu::FindCloset3Mu  -----

  int Tauto3Mu::Get3MuDphi(unsigned target, unsigned obj1, unsigned obj2) {
    if (target == obj1 || obj1 == obj2 || target == obj2)
      return 99999;
    int dPhi1 = deltaPhi(trkMus->at(target).hwPhi(), trkMus->at(obj1).hwPhi());
    int dPhi2 = deltaPhi(trkMus->at(target).hwPhi(), trkMus->at(obj2).hwPhi());
    int dEta1 = deltaEta(trkMus->at(target).hwEta(), trkMus->at(obj1).hwEta());
    int dEta2 = deltaEta(trkMus->at(target).hwEta(), trkMus->at(obj2).hwEta());
    return abs(dPhi1) + abs(dPhi2) + abs(dEta1) + abs(dEta2);
  }
}  // namespace Phase2L1GMT

#endif  // ----- #ifndef PHASE2GMT_TAUTO3MU -----
