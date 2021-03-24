#include <cmath>
#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuTrackTransform.h"

void L1TkEtMissEmuTrackTransform::GenerateLUTs(){
    phi_quadrants = generate_phi_slice_LUTs (L1TkEtMissEmuAlgo::N_quadrants_);
    phi_shift     = generate_phi_slice_LUTs (L1TkEtMissEmuAlgo::N_sectors_);
}

InternalEtWord L1TkEtMissEmuTrackTransform::TransformTrack(TTTrack< Ref_Phase2TrackerDigi_ >& track_ref, float PV){
    InternalEtWord Outword;

    Outword.pV = digitize_Signed <iZ0> (PV, -maxmap["z0"], maxmap["z0"], binmap["z0"]); //Convert vertex 
    //Read in and convert track parameters to integer representation
    Outword.pt   = digitize_Signed <iPt> (track_ref.momentum().perp(),0,maxmap["pt"],binmap["pt"]);
    Outword.eta = digitize_Signed <iEta> (abs(track_ref.momentum().eta()),-maxmap["eta"], maxmap["eta"], binmap["eta"]);

    Outword.chi2rphidof = Chi_to_FWChi <ichi2XY> (track_ref.chi2XY(),chi2Values,binmap["chixy"]);   
    Outword.chi2rzdof   = Chi_to_FWChi <ichi2XY> (track_ref.chi2Z(),chi2Values,binmap["chixy"]);
    Outword.bendChi2  = Chi_to_FWChi <ichi2bend> (track_ref.stubPtConsistency(),chi2bendValues,binmap["chibend"]);
    Outword.nstubs       = CountNStub(track_ref.unpack_hitPattern());

    unsigned int Sector  = track_ref.phiSector();
    Outword.Sector  = Sector;
    // convert to local phi
    Outword.phi = track_ref.phi();
    iPhi localPhi        = FloatPhi_to_iPhi(track_ref.phi(),Sector);
    // Convert to global phi
    Outword.globalPhi   = local_to_global(localPhi,phi_shift[Sector]);

    Outword.z0       = digitize_Signed <iZ0> (track_ref.z0(), -maxmap["z0"], maxmap["z0"], binmap["z0"]);

    return Outword;
}

InternalEtWord L1TkEtMissEmuTrackTransform::TransfromCuts(const edm::ParameterSet& iConfig){
  InternalEtWord cutword;
  cutword.z0          = digitize_Signed <iZ0> ((float)iConfig.getParameter<double>("maxZ0"),-maxmap["z0"], maxmap["z0"], binmap["z0"]);

  cutword.eta        = digitize_Signed <iEta> ((float)iConfig.getParameter<double>("maxEta"), -maxmap["eta"],maxmap["eta"], binmap["eta"]);

  cutword.chi2rphidof = Chi_to_FWChi <ichi2XY> ((float)iConfig.getParameter<double>("chi2rphidofMax"),chi2Values,binmap["chixy"]);
  cutword.chi2rzdof   = Chi_to_FWChi <ichi2XY> ((float)iConfig.getParameter<double>("chi2rzdofMax"),chi2Values,binmap["chixy"]);
  cutword.bendChi2   = Chi_to_FWChi <ichi2bend> ((float)iConfig.getParameter<double>("bendChi2Max"),chi2bendValues,binmap["chibend"]);
  cutword.pt          = digitize_Signed <iPt> ((float)iConfig.getParameter<double>("minPt"),0,maxmap["pt"],binmap["pt"]);
  cutword.nstubs  = (iNstub)iConfig.getParameter<int>("nStubsmin");

  return cutword;

}

iglobPhi L1TkEtMissEmuTrackTransform::local_to_global(iPhi local_phi,iglobPhi sector_shift){ 
      int PhiShift = binmap["globphi"]/2;
      int PhiMin = phi_quadrants.front();
      int PhiMax = phi_quadrants.back();
      int phi_multiplier = L1TkEtMissEmuAlgo::N_phiBits_ - L1TkEtMissEmuAlgo::N_globphiBits_;

      int tempPhi = 0;
      iglobPhi globalPhi = 0;

      tempPhi = (local_phi / pow(2,phi_multiplier)) + sector_shift - PhiShift;
      if (tempPhi < PhiMin) {tempPhi = tempPhi + PhiMax;}
      else if(tempPhi > PhiMax) {tempPhi = tempPhi - PhiMax;}
      else tempPhi = tempPhi;

      globalPhi = iglobPhi (tempPhi);

      return globalPhi;

}

iNstub L1TkEtMissEmuTrackTransform::CountNStub(unsigned int Hitpattern){  
    iNstub Nstub = 0;
    for (int i = (L1TkEtMissEmuAlgo::N_hitpatternbits_-1); i >= 0; i--) {
      int k = Hitpattern >> i;
        if (k & 1)
          Nstub ++;
    }
    return Nstub;
}

iPhi L1TkEtMissEmuTrackTransform::FloatPhi_to_iPhi(float phi,unsigned int sector){   
    float temp_phi = 0.0;
    if (sector < 4 ){ 
      temp_phi = phi - (sector * (2*M_PI)/9); 
    } else if (sector > 5 ){
      temp_phi = phi + ((9-sector) * (2*M_PI)/9); 
    } else if (sector == 4 ){
      if (phi > 0){
        temp_phi = phi - (sector *  (2*M_PI)/9);
      } else{ 
        temp_phi = phi + ((9-sector) * (2*M_PI)/9);
      }
    } else if (sector == 5 ){
      if (phi < 0){
        temp_phi = phi + ((9-sector) * (2*M_PI)/9); 
      } else{
         temp_phi = phi - (sector *  (2*M_PI)/9);
      }
    }
    return L1TkEtMissEmuAlgo::digitize_Signed <iPhi> (temp_phi,-maxmap["phi"],maxmap["phi"],binmap["phi"]);
}

std::vector<iglobPhi> L1TkEtMissEmuTrackTransform::generate_phi_slice_LUTs(unsigned int N) {
    float slice_centre= 0.0;
    std::vector<iglobPhi> phi_LUT;
    for (unsigned int q=0;q<=N;q++){
      phi_LUT.push_back((iglobPhi)(slice_centre / (2*maxmap["phi"] / (binmap["globphi"]-1))));
      slice_centre += 2*M_PI/N;
    }
    return phi_LUT;
}

template <typename chi>
chi L1TkEtMissEmuTrackTransform::Chi_to_FWChi(float Chi2,const float bins[],unsigned int N){  
    chi outputChi = 0;
    for (unsigned int ibin = 0; ibin < N; ++ibin) {
      outputChi = ibin;
      if (Chi2 < bins[ibin])
        break;
    }
    return outputChi;
  }
  