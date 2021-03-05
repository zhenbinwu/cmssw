#ifndef L1Trigger_L1TTrackMatch_L1TkEtMissEmuAlgo_HH
#define L1Trigger_L1TTrackMatch_L1TkEtMissEmuAlgo_HH

#include "DataFormats/L1TrackTrigger/interface/TTTrack_TrackWord.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include <ap_int.h>
#include <math.h>
#include <iostream>
#include <cstdint>
#include <filesystem>
#include <fstream>

using namespace std;
namespace L1TkEtMissEmuAlgo {

  const unsigned int N_chi2XYbits_{4};  //Taken from TTTrack_Word class
  const unsigned int N_chi2bendbits_{3};
  const unsigned int N_hitpatternbits_{7};
  const unsigned int N_ptBits_{15};
  const unsigned int N_etaBits_{16};
  const unsigned int N_z0Bits_{12};
  const unsigned int N_phiBits_{12};

  const unsigned int N_globphiBits_{11};  //Precision of Cos LUT (only spans 0 to pi/2) whereas global phi ( 0 to 2 pi)
  const unsigned int N_globphiExtra_{3};  //Extra bits needed by global phi to span full range
  const unsigned int N_EtExtra_{8};  //Extra room for Et sums
  
  
  typedef ap_uint<N_chi2XYbits_> ichi2XY;
  typedef ap_uint<N_chi2bendbits_> ichi2bend;
  typedef ap_uint<3> iNstub;
  typedef ap_uint<N_phiBits_>  iPhi;
  typedef ap_uint<N_globphiBits_+ N_globphiExtra_>  iglobPhi;
  typedef ap_uint<N_z0Bits_>  iZ0;
  typedef ap_uint<N_ptBits_> iPt;
  typedef ap_int <N_ptBits_+N_EtExtra_> iEt;
  typedef ap_uint<N_etaBits_> iEta;
  
  
  /*
  typedef unsigned int ichi2XY;
  typedef unsigned int ichi2bend;
  typedef unsigned int iNstub;
  typedef unsigned int iPhi;
  typedef unsigned int iglobPhi;
  typedef unsigned int iZ0;
  typedef unsigned int iPt;
  typedef int iEt;
  typedef unsigned int iEta;
  */
  
  const unsigned int N_chi2XYbins_   = 1 << N_chi2XYbits_;
  const unsigned int N_chi2bendbins_ = 1 << N_chi2bendbits_;
  const unsigned int N_ptBins_       = 1 << N_ptBits_;
  const unsigned int N_etaBins_      = 1 << N_etaBits_;
  const unsigned int N_z0Bins_       = 1 << N_z0Bits_;
  const unsigned int N_phiBins_      = 1 << N_phiBits_;
  const unsigned int N_globPhiBins_  = 1 << N_globphiBits_;

  const float max_TWord_Pt_{1024};  //21955               // Arbitrary, define properly!
  const float max_TWord_Phi_{0.69813248};         // relative to the center of the sector
  const float max_TWord_Eta_{2.776472281}; 
  const float max_TWord_Z0_{20.46912512};         

  const unsigned int N_etaregions_{6};
  const unsigned int N_sectors_{9};

  const unsigned int N_quadrants_{4}; //4 quadrants + 1

  const float chi2Values_[N_chi2XYbins_] = {0, 0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 40, 100, 200, 500, 1000, 3000};
  const float chi2bendValues_[N_chi2bendbins_] = {0, 0.5, 1.25, 2, 3, 5, 10, 50};

  const float EtaRegions_[N_etaregions_+1] = {0,0.7,1.0,1.2,1.6,2.0,2.4};
  const float DeltaZ_[N_etaregions_] = {0.4,0.6,0.76,1.0,1.7,2.2};

  const float max_LUT_phi_{M_PI/2};   // Enough symmetry in cos and sin between 0 and pi/2 to get all possible values of cos and sin phi
  const int cosLUT_bins_ = ceil((max_LUT_phi_ * (N_globPhiBins_ -1 ))/(2*max_TWord_Phi_));  //To have same bin spacing between 0 and pi/2 as between original phi granularity

  const string LUTdir_{"LUTs/"};

  struct EtMiss{
    iEt Et;
    iglobPhi Phi;
  };

  //Function to count stubs in hitpattern
  iNstub CountNStub(unsigned int Hitpattern);
  //function to transfrom float chi to int chi bins -> would be part of tttrack, needed for chi cuts
  template <typename chi>
  chi Chi_to_FWChi(float Chi2,const float bins[],unsigned int N){  
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
  T digitize_Signed(float var, float min, float max, unsigned int n_bins) {  
    float temp_var{0.0};
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
  std::vector<iglobPhi> FillCosLUT(unsigned int cosLUT_size);  
  // For generating phi slices for LUTs
  std::vector<iglobPhi> generate_phi_slice_LUTs(unsigned int N);
  
  //Converts local int phi to global int phi
  iglobPhi local_to_global(iPhi local_phi,iglobPhi sector_shift,std::vector<iglobPhi> phi_quadrants);
  // Generate Eta LUT for track to vertex association
  std::vector<iEta> generate_EtaRegions(const float EtaRegions[], unsigned int Nreg );
  // Generate deltaZ bins for track to vertex associaiton
  std::vector<iZ0> generate_DeltaZBins(const float DeltaZ[], unsigned int Nbin );


  template <typename T>
  void writevectorLUTout(vector<T> (&LUT),const string& filename,const string& delimiter){
    int fail = system((string("mkdir -p ") + L1TkEtMissEmuAlgo::LUTdir_).c_str());
    if (fail)
      throw cms::Exception("BadDir") << __FILE__ << " " << __LINE__ << " could not create directory "
                                     << L1TkEtMissEmuAlgo::LUTdir_;
    

    const string fname = L1TkEtMissEmuAlgo::LUTdir_ + filename + ".tab";
    ofstream outstream(fname);
    if (outstream.fail())
      throw cms::Exception("BadFile") << __FILE__ << " " << __LINE__ << " could not create file " << fname;

    outstream << "{" << endl;
    for (unsigned int i = 0; i < LUT.size(); i++) {
      if (i != 0)
        outstream << delimiter << endl;
      outstream << LUT[i];
    }
    outstream << endl << "};" << endl;
    outstream.close();
  }

}  // namespace L1TkEtMissEmuAlgo
#endif