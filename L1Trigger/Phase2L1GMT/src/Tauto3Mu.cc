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

#include "L1Trigger/Phase2L1GMT/interface/Tauto3Mu.h"
#include "TLorentzVector.h"

using namespace Phase2L1GMT;

Tauto3Mu::Tauto3Mu(const edm::ParameterSet &iConfig):
  verbose_(iConfig.getParameter<int>("verbose")),
  hwdRcut_(iConfig.getParameter<int>("pairdRCut")),
  hwdzcut_(iConfig.getParameter<int>("pairdzCut")),
  dodRcut_(iConfig.getParameter<bool>("dopairdRCut")),
  dodzcut_(iConfig.getParameter<bool>("dopairdzCut")),
  dumpForHLS_(iConfig.getParameter<int>("IsodumpForHLS")) {
    if (dumpForHLS_) {
      dumpInput.open("ta3mu_input_data.txt", std::ofstream::out);
      dumpOutput.open("HLS_m3mu_testebench.txt", std::ofstream::out);
    }
    trkMus = new std::vector<l1t::TrackerMuon>();
  }

Tauto3Mu::~Tauto3Mu() {
  if (dumpForHLS_) {
    dumpInput.close();
    dumpOutput.close();
  }
}

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
    trkMus->at(obj0).hwPt() << " "<< trkMus->at(obj0).hwPhi() << " "<< trkMus->at(obj0).hwEta() << " "<< trkMus->at(obj0).hwEta() << " "<<
    trkMus->at(obj1).hwPt() << " "<< trkMus->at(obj1).hwPhi() << " "<< trkMus->at(obj1).hwEta() << " "<< trkMus->at(obj1).hwEta() << " "<<
    trkMus->at(obj2).hwPt() << " "<< trkMus->at(obj2).hwPhi() << " "<< trkMus->at(obj2).hwEta() << " "<< trkMus->at(obj2).hwEta() << " "<<
    expMass<<"  "<< emuMass << endl;
}

void Tauto3Mu::DumpOutputs(int Nevt, unsigned obj0, unsigned obj1, unsigned obj2, float expMass, float emuMass) {
  dumpOutput << Nevt <<" : ";
  DumpOutputs(obj0, obj1, obj2, expMass, emuMass);
}

// ===  FUNCTION  ============================================================
//         Name:  Tauto3Mu::GetTau3Mu
//  Description:
// ===========================================================================
std::vector<l1t::Tau23Mu> Tauto3Mu::GetTau3Mu(const std::vector<l1t::TrackerMuon> &trkMus_) {
  outtau.clear();
  Nevt++;

  if (trkMus_.size() < 3)
    return outtau;

  trkMus->clear();
  for(auto i : trkMus_)
  {
    trkMus->push_back(i);
  }
  Find3MuComb(*trkMus);
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
    if( (nearby3mu[i][0] == nearby3mu[i][1])  || 
        (nearby3mu[i][1] == nearby3mu[i][2])  ||
        (nearby3mu[i][2] == nearby3mu[i][0]) )
      continue;

    float trimass = Get3MuMass(nearby3mu[i][0], nearby3mu[i][1], nearby3mu[i][2]);
    int shiftedmass = int(trimass) >>3;
    ap_uint<8> outmass = ~0;
    if (shiftedmass < outmass )
    {
      outmass = shiftedmass;
    }

    l1t::Tau23Mu temp(seedmu.hwPt(), seedmu.hwEta(), seedmu.hwPhi(), trimass, 
        0, nearby3mu[i][0], nearby3mu[i][1], nearby3mu[i][2]);

    outtau.push_back(temp);

    if (dumpForHLS_) {
      float expMass = DumpInputs(nearby3mu[i][0], nearby3mu[i][1], nearby3mu[i][2]);
      float cvtmass = sqrt(2*trimass) * LSBpt *(1<<PT_MINV_LSB_RED);
      std::cout << " expMass " << expMass <<" cvtMass " << cvtmass << " outmass " << outmass 
        <<" expIntMAss " << cvtmass*cvtmass/0.25 <<" diff " << cvtmass*cvtmass/outmass << std::endl;
      DumpOutputs(Nevt, nearby3mu[i][0], nearby3mu[i][1], nearby3mu[i][2], expMass, outmass);
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

inline int Tauto3Mu::deltaEta(const int eta1, const int eta2) {
  static const int maxbits = (1 << BITSETA) - 1;
  int deta = abs(eta1 - eta2);
  deta &= maxbits;
  return deta;
}
// Ideal the object should carry its own ap types once we finalize
inline int Tauto3Mu::deltaPhi(int phi1, int phi2) {
  static const int maxbits = (1 << BITSPHI) - 1;
  static const int pibits = (1 << (BITSPHI - 1));
  int dphi = abs(phi1 - phi2);
  if (dphi >= pibits)
    dphi = maxbits - dphi;
  return dphi;
}
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
    if (mindR == 99999)
    {
      std::vector<unsigned> tempv={i, i, i};
    }
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

  if (dodRcut_ && ((dPhi1 + dEta1) > hwdRcut_ || (dPhi2 + dEta2) > hwdRcut_ ))
    return 99999;
  if (dodzcut_ && (  abs(trkMus->at(target).hwZ0() - trkMus->at(obj1).hwZ0() ) > hwdzcut_
        || abs(trkMus->at(target).hwZ0() - trkMus->at(obj2).hwZ0() )> hwdzcut_) )
    return 99999;
  return abs(dPhi1) + abs(dPhi2) + abs(dEta1) + abs(dEta2);
}
