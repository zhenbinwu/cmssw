#ifndef L1Trigger_L1TTrackMatch_METTrackTransform_HH
#define L1Trigger_L1TTrackMatch_METTrackTransform_HH

#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuAlgo.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

using namespace L1TkEtMissEmuAlgo;

struct InternalEtWord{
    iZ0 pV;

    iZ0 z0;
    iPt pt;
    iEta eta;
    iglobPhi globalPhi;
    iNstub nstubs;
    ichi2bend bendChi2;
    ichi2XY chi2rphidof;
    ichi2XY chi2rzdof;

    unsigned int Sector;

    float phi;  //Used to debug cos phi LUT
};

class METTrackTransform {
public:
  METTrackTransform() = default;  
  ~METTrackTransform() {};

  void GenerateLUTs();  //Generate internal LUTs needed for track transfrom

  InternalEtWord TransformTrack(TTTrack< Ref_Phase2TrackerDigi_ >& track_ref,float PV);
  InternalEtWord TransfromCuts(const edm::ParameterSet& iConfig);

    //Converts local int phi to global int phi
  iglobPhi local_to_global(iPhi local_phi,iglobPhi sector_shift);

    //Function to count stubs in hitpattern
  iNstub CountNStub(unsigned int Hitpattern);

  //Function to take float phi to local integer phi -> would be part of tttrack
  iPhi FloatPhi_to_iPhi(float phi,unsigned int sector);

  std::vector<iglobPhi> generate_phi_slice_LUTs(unsigned int N);

  //function to transfrom float chi to int chi bins -> would be part of tttrack, needed for chi cuts
  template <typename chi>
  chi Chi_to_FWChi(float Chi2,const float bins[],unsigned int N);

  std::vector<iglobPhi> get_phi_quad() const {return phi_quadrants;}
  std::vector<iglobPhi> get_phi_shift() const {return phi_shift;}

  std::map<std::string, unsigned int> get_binmap() const {return binmap;}
  std::map<std::string, float> get_maxmap() const {return maxmap;}


private:

  std::map<std::string, unsigned int> binmap { {"chixy",N_chi2XYbins_}, 
                                               {"chibend",N_chi2bendbins_},
                                               {"pt",N_ptBins_},
                                               {"eta",N_etaBins_} ,
                                               {"z0",N_z0Bins_},
                                               {"phi",N_phiBins_},
                                               {"globphi",N_globPhiBins_}};

  std::map<std::string, float> maxmap { {"pt",max_TWord_Pt_},
                                        {"eta",max_TWord_Eta_},
                                        {"z0",max_TWord_Z0_},
                                        {"phi",max_TWord_Phi_}};


  std::vector<iglobPhi> phi_quadrants;
  std::vector<iglobPhi> phi_shift;

  const float chi2Values[N_chi2XYbins_] = {0, 0.25, 0.5, 1, 2, 3, 5, 7, 10, 20, 40, 100, 200, 500, 1000, 3000};
  const float chi2bendValues[N_chi2bendbins_] = {0, 0.5, 1.25, 2, 3, 5, 10, 50};




 
};

#endif