#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuAlgo.h"

namespace L1TkEtMissEmuAlgo{
  //Function to count stubs in hitpattern
  iNstub CountNStub(unsigned int Hitpattern){  
    iNstub Nstub = 0;
    for (int i = (N_hitpatternbits_-1); i >= 0; i--) {
      int k = Hitpattern >> i;
        if (k & 1)
          Nstub ++;
    }
    return Nstub;
  }

  //Function to take float phi to local integer phi -> would be part of tttrack
  iPhi FloatPhi_to_iPhi(float phi,unsigned int sector){   
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
    return digitize_Signed <iPhi> (temp_phi,-max_TWord_Phi_,max_TWord_Phi_,N_phiBins_);
  }


  std::vector<iglobPhi> FillCosLUT(unsigned int cosLUT_size){  //Fill cosine LUT with integer values
    float phi = 0;
    std::vector<iglobPhi> cosLUT;
    for (unsigned int LUT_idx = 0; LUT_idx <= cosLUT_size; LUT_idx++){
      cosLUT.push_back(digitize_Signed <iglobPhi> (cos(phi),0,1,N_globPhiBins_));
      phi += (2*max_TWord_Phi_)/(N_globPhiBins_ -1 );
     }
    return cosLUT;
  }

  std::vector<iglobPhi> generate_phi_slice_LUTs(unsigned int N) {
    float slice_centre= 0.0;
    std::vector<iglobPhi> phi_LUT;
    for (unsigned int q=0;q<=N;q++){
      phi_LUT.push_back((iglobPhi)(slice_centre / (2*max_TWord_Phi_ / (N_globPhiBins_-1))));
      slice_centre += 2*M_PI/N;
    }
    return phi_LUT;
  }

  //Converts local int phi to global int phi
  iglobPhi local_to_global(iPhi local_phi,iglobPhi sector_shift,std::vector<iglobPhi> phi_quadrants){ 
      int PhiShift = N_globPhiBins_/2;
      int PhiMin = phi_quadrants.front();
      int PhiMax = phi_quadrants.back();
      int phi_multiplier = N_phiBits_ - N_globphiBits_;

      int tempPhi = 0;
      iglobPhi globalPhi = 0;

      tempPhi = (local_phi / pow(2,phi_multiplier)) + sector_shift - PhiShift;
      if (tempPhi < PhiMin) {tempPhi = tempPhi + PhiMax;}
      else if(tempPhi > PhiMax) {tempPhi = tempPhi - PhiMax;}
      else tempPhi = tempPhi;

      globalPhi = iglobPhi (tempPhi);

      return globalPhi;

  }

    // Generate Eta LUT for track to vertex association
  std::vector<iEta> generate_EtaRegions(const float EtaRegions[],unsigned int Nreg) {
    std::vector<iEta> LUT;
    for (unsigned int q=0;q<Nreg;q++){
      LUT.push_back(digitize_Signed <iEta> (EtaRegions[q], -max_TWord_Eta_, max_TWord_Eta_, N_etaBins_));
    }
    return LUT;
    
  }

  std::vector<iZ0> generate_DeltaZBins(const float DeltaZ[],unsigned int Nbin) {
    std::vector<iZ0> LUT;
    for (unsigned int q=0;q<Nbin;q++){
      LUT.push_back(digitize_Signed <iZ0> (DeltaZ[q], 0, max_TWord_Z0_, N_z0Bins_/2));  //Only half range otherwise deltaZ values midway through integer representation
    }
    return LUT;
    
  }

}