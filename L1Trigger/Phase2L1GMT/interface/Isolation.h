// ===========================================================================
// 
//       Filename:  Isolation.h
// 
//    Description:  
// 
//        Version:  1.0
//        Created:  02/23/2021 01:16:43 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Zhenbin Wu (benwu), zhenbin.wu@gmail.com
//        Company:  UIC, CMS@LPC, CDF@FNAL
// 
// ===========================================================================

#ifndef  PHASE2GMT_ISOLATION
#define  PHASE2GMT_ISOLATION

#include "DataFormats/L1TMuonPhase2/interface/TrackerMuon.h"
#include "L1Trigger/Phase2L1GMT/interface/ConvertedTTTrack.h"

#include <fstream>

namespace Phase2L1GMT
{
  class Isolation
  {
    public:
      Isolation(const edm::ParameterSet& iConfig);
      ~Isolation();
      Isolation( const Isolation &cpy );

      unsigned compute_trk_iso( l1t::TrackerMuon &in_mu, ConvertedTTTrack &in_trk);

      void isolation_allmu_alltrk( std::vector<l1t::TrackerMuon> &trkMus, 
          std::vector<ConvertedTTTrack> &convertedTracks);

    private:
      void DumpOutputs( std::vector<l1t::TrackerMuon> &trkMus);
      void DumpInputs( std::vector<l1t::TrackerMuon> &trkMus, std::vector<ConvertedTTTrack> &convertedTracks);

      int deltaEta(const int eta1,const int eta2);
      int deltaZ0(const int Z01,const int Z02);
      int deltaPhi(int phi1, int phi2);


      // Current setting in Luca's code
      const static int BITMUPT          = 14;
      const static int BITMUPHI         = 13;
      const static int BITMUETA         = 13;
      const static int BITMUZ0          = 12;
      const static int BITMUD0          = 13;
      const static int BITMUCHARGE      = 1;
      const static int BITMUQUALITY     = 4;
      const static int BITMUISO         = 2;

      const double lsb_pt  = 512./pow(2, BITMUPT);
      const double lsb_eta = 2.*M_PI/pow(2, BITMUETA);
      const double lsb_phi = 2.*M_PI/pow(2, BITMUPHI);
      const double lsb_z0  = 60./pow(2, BITMUZ0);
      const static int c_iso_dangle_max = 260;//@ <  260 x 2pi/2^13 = 0.2 rad
      const static int c_iso_dz_max     = 68; //@ <  68  x 60/2^12  = 1   cm
      const static int c_iso_pt_min     = 96; //@ >= , 96  x 1/2^5 = 3 GeV
      // Assuming 4 bits for Muon isolation
      int iso_threshold_1;
      int iso_threshold_2;
      int iso_threshold_3;
      int iso_threshold_4;
      bool verbose_;
      bool dumpForHLS_;
      std::ofstream dumpInput;
      std::ofstream dumpOutput;

  };

  Isolation::Isolation(const edm::ParameterSet& iConfig):
    iso_threshold_1(iConfig.getParameter<int>("IsoThreshold1")),
    iso_threshold_2(iConfig.getParameter<int>("IsoThreshold2")),
    iso_threshold_3(iConfig.getParameter<int>("IsoThreshold3")),
    iso_threshold_4(iConfig.getParameter<int>("IsoThreshold4")),
    verbose_(iConfig.getParameter<int>("verbose")),
    dumpForHLS_(iConfig.getParameter<int>("IsodumpForHLS"))
  {
    dumpForHLS_ = true;
    if (dumpForHLS_)
    {
      dumpInput.open("Isolation_Mu_Track_infolist.txt", std::ofstream::out);
      dumpOutput.open("Isolation_Mu_Isolation.txt", std::ofstream::out);
    }
  }

  Isolation::~Isolation()
  {

    if (dumpForHLS_)
    {
      dumpInput.close();
      dumpOutput.close();
    }
  }

  Isolation::Isolation( const Isolation &cpy )
  {
  }


  void Isolation::DumpInputs( std::vector<l1t::TrackerMuon> &trkMus, std::vector<ConvertedTTTrack>  &convertedTracks)
  {
    static int nevti =0;
    nevti++;
    int totalsize =0;
    int exptotal = 12 + 18 * 100; // N_Muon + N_TRK_LINKS * NTRKperlinks
    for (unsigned int i = 0; i < trkMus.size(); ++i)
    {
      dumpInput <<" " << nevti <<" 0 " << i <<" "<<trkMus.at(i).hwPt() * lsb_pt<<" " <<trkMus.at(i).hwEta() * lsb_eta
        <<" " << trkMus.at(i).hwPhi() * lsb_phi <<" " << trkMus.at(i).hwZ0() * lsb_z0<< std::endl;
      totalsize++;
    }
    for (unsigned int i = 0; i < convertedTracks.size(); ++i)
    {
      dumpInput <<" " << nevti <<" 1 "<< i <<" "<<convertedTracks.at(i).pt() * lsb_pt<<" " <<convertedTracks.at(i).eta() * lsb_eta
        <<" " << convertedTracks.at(i).phi() * lsb_phi<<" " << convertedTracks.at(i).z0() *lsb_z0<< std::endl;
      totalsize++;
    }
    int ntrks = convertedTracks.size();
    // Pat the remaining 
    while (totalsize  < exptotal)
    {
      dumpInput<<" " << nevti << " 1 " << ntrks++ <<" " <<0 <<" " <<0 <<" " <<0 <<" " <<0  <<" " <<0 <<" " <<0 <<std::endl;
      totalsize++;
    }
  }

  void Isolation::DumpOutputs( std::vector<l1t::TrackerMuon> &trkMus)
  {
    static int nevto =0;
    for (unsigned int i = 0; i < trkMus.size(); ++i)
    {
      auto mu = trkMus.at(i);
      if (mu.hwPt() != 0)
      {
        double convertphi = mu.hwPhi() *lsb_phi;
        if (convertphi > M_PI)
        {
          convertphi -= 2 * M_PI;
        }
        dumpOutput << nevto << " " << i << " " << mu.hwPt() *lsb_pt << " " << mu.hwEta() *lsb_eta<< " " 
          << convertphi << " " << mu.hwZ0() *lsb_z0  << " " << mu.hwIso() << endl;
      }
    }
    nevto++;
  }

  void Isolation::isolation_allmu_alltrk( std::vector<l1t::TrackerMuon> &trkMus, std::vector<ConvertedTTTrack>  &convertedTracks)
  {

    if(dumpForHLS_)
    {
      DumpInputs(trkMus, convertedTracks);
    }

    for(auto &mu : trkMus)
    {
      int accum = 0;
      int iso_ = 0;
      for(auto t : convertedTracks)
      {
        accum += compute_trk_iso(mu, t);
      }

      // Only 8 bit for accumation? 
      ap_ufixed<8, 8, AP_TRN, AP_SAT> iso_accum_t(accum);
      accum = iso_accum_t.to_int();

      iso_ |= (accum < iso_threshold_1) << 0;
      iso_ |= (accum < iso_threshold_2) << 1;
      iso_ |= (accum < iso_threshold_3) << 2;
      iso_ |= (accum < iso_threshold_4) << 3;
      mu.setHwIso(iso_);
    }

    if(dumpForHLS_)
    {
      DumpOutputs(trkMus);
    }

  }

  int Isolation::deltaEta(const int eta1,const int eta2) {
    ap_int<BITMUETA> ap_eta1(eta1);
    ap_int<BITMUETA> ap_eta2(eta2);
    ap_int<BITMUETA> dEta = ap_eta1 - ap_eta2;

    if (dEta<0)
      return (-1*dEta).to_int();
    else
      return dEta.to_int();
  }
  int Isolation::deltaZ0(const int Z01,const int Z02) {
    ap_int<BITMUZ0> ap_Z01(Z01);
    ap_int<BITMUZ0> ap_Z02(Z02);
    ap_int<BITMUZ0> dZ0 = ap_Z01 - ap_Z02;

    if (dZ0<0)
      return (-1*dZ0).to_int();
    else
      return dZ0.to_int();
  }

  // Bad idea here, but as a temparory workout
  // Ideal the object should carry its own ap types once we finalize
  int Isolation::deltaPhi(int phi1, int phi2) {
    ap_int<BITMUPHI> ap_phi1(phi1);
    ap_int<BITMUPHI> ap_phi2(phi2);
    ap_int<BITMUPHI> dPhi = phi1 - phi2;
    ap_uint<BITMUPHI-1> dPhi_unsigned;
    if (dPhi<0)
      dPhi_unsigned = -dPhi;
    else
      dPhi_unsigned = dPhi;
    return dPhi_unsigned.to_int();
  }

  unsigned Isolation::compute_trk_iso(l1t::TrackerMuon  &in_mu, ConvertedTTTrack &in_trk)
  {
    int dphi = deltaPhi(in_mu.hwPhi(),  in_trk.phi());
    int deta = deltaEta(in_mu.hwEta(), in_trk.eta());
    int dz0 = deltaZ0(in_mu.hwZ0(), in_trk.z0());

    // Need to convert, but still figure out the input format. Ignoring it for now
    //bool pass_deta        = (deta      < deta_t(c_iso_dangle_max) ? true : false);
    //bool pass_dphi        = (dphi      < dphi_t(c_iso_dangle_max) ? true : false);
    //bool pass_dz0         = (dz0       < dz0_t(c_iso_dz_max)      ? true : false);
    //bool pass_trkpt       = (in_trk.pt() >= hw_pt_t(c_iso_pt_min) ? true : false);
    bool pass_deta        = (deta      < c_iso_dangle_max ? true : false);
    bool pass_dphi        = (dphi      < c_iso_dangle_max ? true : false);
    bool pass_dz0         = (dz0       < c_iso_dz_max     ? true : false);
    bool pass_trkpt       = (in_trk.pt() >= c_iso_pt_min  ? true : false);
    bool pass_ovrl        = (deta > 0 || dphi > 0         ? true : false);

    // match conditions
    if ( pass_deta  &&
        pass_dphi  &&
        pass_dz0   &&
        pass_trkpt &&
        pass_ovrl
       ){
      return in_trk.pt();
    }
    else{
      return 0;
    }
  }
}
#endif   // ----- #ifndef PHASE2GMT_ISOLATION -----
