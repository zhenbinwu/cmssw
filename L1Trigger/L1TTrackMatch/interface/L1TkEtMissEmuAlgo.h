#ifndef L1Trigger_L1TTrackMatch_L1TkEtMissEmuAlgo_HH
#define L1Trigger_L1TTrackMatch_L1TkEtMissEmuAlgo_HH

#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include <ap_int.h>
#include <math.h>


namespace L1TkEtMissEmuAlgo {

  const unsigned int N_chi2XYbits     = 4;
  const unsigned int N_chi2bendbits   = 3;
  const unsigned int N_hitpatternbits = 7;
  const unsigned int N_ptBits         = 15;
  const unsigned int N_etaBits        = 16;
  const unsigned int N_z0Bits         = 12;
  const unsigned int N_phiBits        = 11;


  const unsigned int N_globphiBits    = 8;  //Precision of Cos LUT (only spans 0 to pi/2) whereas global phi ( 0 to 2 pi)
  const unsigned int N_globphiExtra   = 3;  //Extra bits needed by global phi to span full range
  const unsigned int N_EtExtra        = 8;  //Extra room for Et sums
  
  /*
  typedef ap_uint<N_chi2XYbits> ichi2XY;
  typedef ap_uint<N_chi2bendbits> ichi2bend;
  typedef ap_uint<3> iNstub;
  typedef ap_uint<N_phiBits>  iPhi;
  typedef ap_uint<N_globphiBits+N_globphiExtra>  iglobPhi;
  typedef ap_uint<N_z0Bits>  iZ0;
  typedef ap_uint<N_ptBits> iPt;
  typedef ap_int<N_ptBits+N_EtExtra> iEt;
  typedef ap_uint<N_etaBits> iEta;
  */

  typedef unsigned int ichi2XY;
  typedef unsigned int ichi2bend;
  typedef unsigned int iNstub;
  typedef unsigned int iPhi;
  typedef unsigned int iglobPhi;
  typedef unsigned int iZ0;
  typedef unsigned int iPt;
  typedef int iEt;
  typedef unsigned int iEta;

  
  const int N_chi2XYbins =   1 << N_chi2XYbits;
  const int N_chi2bendbins = 1 << N_chi2bendbits;
  const int N_ptBins =       1 << N_ptBits;
  const int N_etaBins =      1 << N_etaBits;
  const int N_z0Bins =       1 << N_z0Bits;
  const int N_phiBins =      1 << N_phiBits;
  const int N_globPhiBins =  1 << N_globphiBits;

  const float maxPt = 1024;                // Arbitrary, define properly!
  const float maxPhi = 0.69813248;         // relative to the center of the sector
  const float maxEta = 2.644120761058629;  // Arbitrary, define properly!
  const float maxZ0 = 20.46912512;         

  const int N_etaregions = 6;
  
  const float chi2Values[N_chi2XYbins] = {0, 0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 40, 100, 200, 500, 1000, 3000};
  const float chi2bendValues[N_chi2bendbins] = {0, 0.5, 1.25, 2, 3, 5, 10, 50};

  const float EtaRegions[N_etaregions+1] = {0,0.7,1.0,1.2,1.6,2.0,2.4};
  const float DeltaZ[N_etaregions] = {0.4,0.6,0.76,1.0,1.7,2.2};



  const float max_LUT_phi =  M_PI/2;   // Enough symmetry in cos and sin between 0 and pi/2 to get all possible values of cos and sin phi
  const int LUT_bins = (max_LUT_phi * (N_globPhiBins - 1))/(2*maxPhi);  //To have same bin spacing between 0 and pi/2 as between original phi granularity

  const int cordic_steps = 4;  //Steps for cordic, higher means more accuracy
  const iPhi phiScale = N_globPhiBins;  //Phi scale for cordic, degrees in radians don't work in integer, too large and overflow phi

  //Function to count stubs in hitpattern
  iNstub CountNStub(unsigned int Hitpattern);
  //function to transfrom float chi to int chi bins -> would be part of tttrack, needed for chi cuts
  template <typename chi,int N>
  chi Chi_to_FWChi(float Chi2,const float bins[N]){  
    chi outputChi = 0;
    for (unsigned int ibin = 0; ibin < N; ++ibin) {
      outputChi = ibin;
      if (Chi2 < bins[ibin])
        break;
    }
    return outputChi;
  }

  //Function to transform float variables to fixed type ints -> would be part of tttrack needed for cuts
  template <typename T>
  T digitize_Signed(float var, float min, float max, int n_bins) {  
    float temp_var = 0.;
    T seg = 0;
    
    if (var > max ) {temp_var = max;}  else {temp_var = var;}
    if (var < min) {temp_var = min;}  else {temp_var = temp_var;}

    temp_var = (temp_var - min) / ((max-min)/(n_bins-1));
    seg = T (temp_var); 
    return seg;
  }

  //Function to take float phi to local integer phi -> would be part of tttrack
  iPhi FloatPhi_to_iPhi(float phi,unsigned int sector);  
  //Fill cosine LUT with integer values
  void FillCosLUT(iglobPhi cosLUT[LUT_bins],float max_LUT_phi);  
  // For generating phi slices for LUTs
  std::vector<iglobPhi> generate_phi_slice_LUTs(unsigned int Num_entries);
  //Converts local int phi to global int phi
  iglobPhi local_to_global(iPhi local_phi,iglobPhi sector_shift,std::vector<iglobPhi> phi_quadrants);
  // Generate Eta LUT for track to vertex association
  std::vector<iEta> generate_EtaRegions(const float EtaRegions[N_etaregions+1]);
  // Generate deltaZ bins for track to vertex associaiton
  std::vector<iZ0> generate_DeltaZBins(const float DeltaZ[N_etaregions]);

  struct EtMiss{
    iEt Et;
    iglobPhi Phi;
  };

}  // namespace L1TkEtMissEmuAlgo
#endif