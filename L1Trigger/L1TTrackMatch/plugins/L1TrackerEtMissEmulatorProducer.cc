// Original Author:  Emmanuelle Perez,40 1-A28,+41227671915,
//         Created:  Tue Nov 12 17:03:19 CET 2013
//Modified by Emily MacDonald, 30 Nov 2018

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TCorrelator/interface/TkEtMiss.h"
#include "DataFormats/L1TCorrelator/interface/TkEtMissFwd.h"
#include "DataFormats/L1TCorrelator/interface/TkPrimaryVertex.h"
#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuAlgo.h"
#include "L1Trigger/L1TTrackMatch/interface/Cordic.h"

#include <numeric>
#include <ap_int.h>

using namespace L1TkEtMissEmuAlgo;
using namespace l1t;

class L1TrackerEtMissEmulatorProducer : public edm::stream::EDProducer<> {
public:
  typedef TTTrack< Ref_Phase2TrackerDigi_ >  L1TTTrackType;
  typedef std::vector< L1TTTrackType > L1TTTrackCollectionType;

  explicit L1TrackerEtMissEmulatorProducer(const edm::ParameterSet&);
  ~L1TrackerEtMissEmulatorProducer();

private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
  iZ0       maxZ0_;	 // in HWU
  iZ0       minZ0_;	 
  iEta      maxEta_; 
  iEta      minEta_; 
  ichi2XY   chi2rphidofMax_; 
  ichi2XY   chi2rzdofMax_; 
  ichi2bend bendChi2Max_;   
  iPt       minPt_; 
  iNstub    nStubsmin_;

  std::vector<iglobPhi> cosLUT ;      // Cos LUT array
  std::vector<iEta> EtaRegionsLUT;    //Various precomputed LUTs
  std::vector<iZ0>  DeltaZLUT;        
  std::vector<iglobPhi> phi_quadrants;
  std::vector<iglobPhi> phi_shift;   

  int CordicSteps;
  int debug;
  int phiScale_;
  bool cordicdebug = false;
  bool writeLUTs   = false;

  std::string L1MetCollectionName;

  std::map<std::string, unsigned int> binmap { {"chixy",N_chi2XYbins_}, 
                                               {"chibend",N_chi2bendbins_},
                                               {"pt",N_ptBins_},
                                               {"eta",N_etaBins_} ,
                                               {"z0",N_z0Bins_},
                                               {"phi",N_phiBins_},
                                               {"globphi",N_globPhiBins_}};

  std::map<std::string, float> maxvalues     { { "pt",  max_TWord_Pt_},
                                               { "phi", max_TWord_Phi_},
                                               { "eta", max_TWord_Eta_},
                                               { "z0",  max_TWord_Z0_}}; 

  const edm::EDGetTokenT< TkPrimaryVertexCollection > pvToken_;
  const edm::EDGetTokenT<std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > > trackToken_;
};

//constructor//
L1TrackerEtMissEmulatorProducer::L1TrackerEtMissEmulatorProducer(const edm::ParameterSet& iConfig) :
pvToken_(consumes<TkPrimaryVertexCollection>(iConfig.getParameter<edm::InputTag>("L1VertexInputTag"))),
trackToken_(consumes< std::vector<TTTrack< Ref_Phase2TrackerDigi_> > > (iConfig.getParameter<edm::InputTag>("L1TrackInputTag")))
{
  // Input parameter cuts and convert to correct integer representations
  maxZ0_          = digitize_Signed <iZ0> ((float)iConfig.getParameter<double>("maxZ0"),-maxvalues["z0"], maxvalues["z0"], binmap["z0"]);
  minZ0_          = digitize_Signed <iZ0> (-(float)iConfig.getParameter<double>("maxZ0"),-maxvalues["z0"], maxvalues["z0"], binmap["z0"]);
  maxEta_         = digitize_Signed <iEta> ((float)iConfig.getParameter<double>("maxEta"), -maxvalues["eta"], maxvalues["eta"], binmap["eta"]);
  minEta_         = digitize_Signed <iEta> (-(float)iConfig.getParameter<double>("maxEta"), -maxvalues["eta"], maxvalues["eta"], binmap["eta"]);
  chi2rphidofMax_ = Chi_to_FWChi <ichi2XY> ((float)iConfig.getParameter<double>("chi2rphidofMax"),chi2Values_,binmap["chixy"]);
  chi2rzdofMax_   = Chi_to_FWChi <ichi2XY> ((float)iConfig.getParameter<double>("chi2rzdofMax"),chi2Values_,binmap["chixy"]);
  bendChi2Max_    = Chi_to_FWChi <ichi2bend> ((float)iConfig.getParameter<double>("bendChi2Max"),chi2bendValues_,binmap["chibend"]);
  minPt_          = digitize_Signed <iPt> ((float)iConfig.getParameter<double>("minPt"),0,maxvalues["pt"],binmap["pt"]);
  nStubsmin_      = (iNstub)iConfig.getParameter<int>("nStubsmin");
  
  CordicSteps     = (int)iConfig.getParameter<int>("nCordicSteps");
  debug           = (int)iConfig.getParameter<int>("Debug");
  writeLUTs       = (bool)iConfig.getParameter<bool>("WriteLUTs");
  phiScale_       = (int)iConfig.getParameter<int>("phiScale");

  L1MetCollectionName = (std::string)iConfig.getParameter<std::string>("L1MetCollectionName");

  if (debug == 5){cordicdebug=true;}

  // Compute LUTs
  cosLUT        = FillCosLUT(cosLUT_bins_);
  phi_quadrants = generate_phi_slice_LUTs (N_quadrants_);
  phi_shift     = generate_phi_slice_LUTs (N_sectors_);
  EtaRegionsLUT = generate_EtaRegions(EtaRegions_,N_etaregions_+1);
  DeltaZLUT     = generate_DeltaZBins(DeltaZ_,N_etaregions_);

  //Write Out LUTs
  if (writeLUTs){
    writevectorLUTout<iglobPhi> (cosLUT,"cos",",");
    writevectorLUTout<iglobPhi> (phi_quadrants,"phiquadrants",",");
    writevectorLUTout<iglobPhi> (phi_shift,"phishift",",");
    writevectorLUTout<iEta>     (EtaRegionsLUT,"etaregions",",");
    writevectorLUTout<iZ0>      (DeltaZLUT,"dzbins",",");
    std::vector<unsigned int> cutlut = {(unsigned int)minPt_,(unsigned int)maxZ0_,
                                        (unsigned int)maxEta_,(unsigned int)chi2rphidofMax_,
                                        (unsigned int)chi2rzdofMax_,(unsigned int)bendChi2Max_,
                                        (unsigned int)nStubsmin_};
    writevectorLUTout<unsigned int>      (cutlut,"cuts",",");
    }

  if (debug == 1){
    edm::LogInfo("L1TkEtMissEmulator") << "=================================================\n";
    edm::LogInfo("L1TkEtMissEmulator") << "Phi Quadrants: ";
    for (auto i: phi_quadrants){ edm::LogInfo("L1TkEtMissEmulator") << i << "|";}
    edm::LogInfo("L1TkEtMissEmulator") << "\n";

    edm::LogInfo("L1TkEtMissEmulator") << "Phi Sector Shift: ";
    for (auto i: phi_shift){ edm::LogInfo("L1TkEtMissEmulator") << i << "|";}
    edm::LogInfo("L1TkEtMissEmulator") << "\n";
    }

  
  produces<TkEtMissCollection>(L1MetCollectionName);
}

L1TrackerEtMissEmulatorProducer::~L1TrackerEtMissEmulatorProducer() { }

void L1TrackerEtMissEmulatorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::unique_ptr<TkEtMissCollection> METCollection(new TkEtMissCollection);

  edm::Handle<TkPrimaryVertexCollection> L1VertexHandle;
  iEvent.getByToken(pvToken_,L1VertexHandle);

  edm::Handle<L1TTTrackCollectionType> L1TTTrackHandle;
  iEvent.getByToken(trackToken_, L1TTTrackHandle);

  Cordic cordic_sqrt(phiScale_,N_ptBits_, CordicSteps,cordicdebug,writeLUTs);


  if( !L1VertexHandle.isValid() ) {
    LogError("L1TrackerEtMissEmulatorProducer")
    << "\nWarning: TkPrimaryVertexCollection not found in the event. Exit\n";
    return;
  }

  if( !L1TTTrackHandle.isValid() ) {
    LogError("L1TrackerEtMissEmulatorProducer")
    << "\nWarning: L1TTTrackCollection not found in the event. Exit\n";
    return;
  }

  // Initialize sector sums, need 0 initialization in case a sector has no tracks
  iEt sumPx[N_sectors_] = {  0,0,0,0,0,0,0,0,0 };
  iEt sumPy[N_sectors_] = {  0,0,0,0,0,0,0,0,0 };

  // Track counters
  int num_tracks{0};
  int num_assoc_tracks{0};
  int num_quality_tracks{0};

  iZ0 zVTX = digitize_Signed <iZ0> (L1VertexHandle->begin()->zvertex(), -maxvalues["z0"], maxvalues["z0"], binmap["z0"]); //Convert vertex TODO use emulator

  for (const auto& track : *L1TTTrackHandle) {
    num_tracks++;
    L1TTTrackType& track_ref = const_cast< L1TTTrackType&>(track); //Non Const member functions in TTTrack_TrackWord 

    //Read in and convert track parameters to integer representation
    iPt pt   = digitize_Signed <iPt> (track_ref.momentum().perp(),0,maxvalues["pt"],binmap["pt"]);
    iEta eta = digitize_Signed <iEta> (abs(track_ref.momentum().eta()),-maxvalues["eta"], maxvalues["eta"], binmap["eta"]);

    ichi2XY chi2rphidof = Chi_to_FWChi <ichi2XY> (track_ref.chi2XY(),chi2Values_,binmap["chixy"]);   
    ichi2XY chi2rzdof   = Chi_to_FWChi <ichi2XY> (track_ref.chi2Z(),chi2Values_,binmap["chixy"]); //TTTrack_TrackWord has no chi2Z getter function 
    ichi2bend bendChi2  = Chi_to_FWChi <ichi2bend> (track_ref.stubPtConsistency(),chi2bendValues_,binmap["chibend"]);
    iNstub nstubs       = CountNStub(track_ref.unpack_hitPattern());

    unsigned int Sector  = track_ref.phiSector();
    // convert to local phi, functionality not in TTTrack_TrackWord word
    iPhi localPhi        = FloatPhi_to_iPhi(track_ref.phi(),Sector);
    // Convert to global phi
    iglobPhi globalPhi   = local_to_global(localPhi,phi_shift[Sector],phi_quadrants);

    iZ0  z0       = digitize_Signed <iZ0> (track_ref.z0(), -maxvalues["z0"], maxvalues["z0"], binmap["z0"]);
    iZ0  deltaZ0  = 0;

    // Parameter cuts
    if (pt <= minPt_) continue;  
    if (z0 > maxZ0_) continue;
    if (z0 < minZ0_) continue;
    if (eta > maxEta_) continue;
    if (eta < minEta_) continue;
    if (chi2rphidof >= chi2rphidofMax_) continue;
    if (chi2rzdof >= chi2rzdofMax_) continue;
    if (bendChi2 >= bendChi2Max_) continue;
    if (nstubs < nStubsmin_) continue;


    num_quality_tracks++;
    // Temporary deltaZ to facilitate use of unsigned int for iZ0
    int temp = z0 - zVTX;
    iZ0 z_diff = abs(temp);
     
    // Track to vertex association
    /*
    if      ( eta>=EtaRegionsLUT[0] &&  eta< EtaRegionsLUT[1])  deltaZ0 = DeltaZLUT[0];
    else if ( eta>=EtaRegionsLUT[1] &&  eta< EtaRegionsLUT[2])  deltaZ0 = DeltaZLUT[1];
    else if ( eta>=EtaRegionsLUT[2] &&  eta< EtaRegionsLUT[3])  deltaZ0 = DeltaZLUT[2];
    else if ( eta>=EtaRegionsLUT[3] &&  eta< EtaRegionsLUT[4])  deltaZ0 = DeltaZLUT[3];
    else if ( eta>=EtaRegionsLUT[4] &&  eta< EtaRegionsLUT[5])  deltaZ0 = DeltaZLUT[4];
    else if ( eta>=EtaRegionsLUT[5] &&  eta<=EtaRegionsLUT[6])  deltaZ0 = DeltaZLUT[5];
    */

    for (unsigned int reg=0;reg < N_etaregions_; reg++){
      if ( eta>=EtaRegionsLUT[reg] &&  eta< EtaRegionsLUT[reg+1]) {
        deltaZ0 = DeltaZLUT[reg];
        break;
      } 
    }

    if (debug == 3){
      std::cout << "========================Z0 Debug=================================\n";
      std::cout << "z0: " << z0 << " z vertex: " << zVTX << '\n';
      std::cout << "Zdiff: " << z_diff << " dZ bin: " << deltaZ0 << '\n';
      std::cout << "eta: " << eta << '\n';
    }

    if ( z_diff <= deltaZ0 ) {
      num_assoc_tracks++;

      if (debug == 2){
        std::cout << "========================Phi Debug=================================\n";
        std::cout << "pt: " << pt << "\n";
        std::cout << "Sector Phi: " << localPhi << " Global Phi: " << globalPhi << "\n";
        std::cout << "Int Phi: "    << globalPhi << " Float Phi: " << track_ref.phi() 
                                               << " Float Cos(Phi): " << cos(track_ref.phi()) 
                                               << " Float Sin(Phi): " << sin(track_ref.phi()) << "\n";
      }
      
      // Split tracks in phi quadrants and access cosLUT, backwards iteration through cosLUT gives sin
      // Sum sector Et -ve when cos or sin phi are -ve
      if (globalPhi >= phi_quadrants[0] && globalPhi < phi_quadrants[1]){
        sumPx[Sector] = sumPx[Sector] + (pt*cosLUT[globalPhi] )/ binmap["globphi"]; 
        sumPy[Sector] = sumPy[Sector] + (pt*cosLUT[phi_quadrants[1] - 1 - globalPhi])/ binmap["globphi"];

        if (debug == 2){
          std::cout << "Sector: " << Sector << " Quadrant: " << 1 << "\n";
          std::cout << "Int Phi: " << globalPhi 
                                                 << " Int Cos(Phi): " << cosLUT[globalPhi] 
                                                 << " Int Sin(Phi): " << cosLUT[phi_quadrants[1]- 1 - globalPhi] << "\n";
          std::cout << "Int Phi: " << (float)globalPhi/binmap["globphi"] 
                                                 << " Int Cos(Phi): " << (float)cosLUT[globalPhi]/binmap["globphi"] 
                                                 << " Int Sin(Phi): " << (float)cosLUT[phi_quadrants[1] - 1- globalPhi]/binmap["globphi"] << "\n";  
        }   
      }
      else if (globalPhi >= phi_quadrants[1] && globalPhi < phi_quadrants[2]){
        sumPx[Sector] = sumPx[Sector] - (pt*cosLUT[phi_quadrants[2] - 1 - globalPhi] )/ binmap["globphi"]; 
        sumPy[Sector] = sumPy[Sector] + (pt*cosLUT[globalPhi - phi_quadrants[1]] )/ binmap["globphi"];  

        if (debug == 2){
          std::cout << "Sector: " << Sector << " Quadrant: " << 2 << "\n";
          std::cout << "Int Phi: " << globalPhi 
                                                 << " Int Cos(Phi): -" << cosLUT[phi_quadrants[2] - globalPhi] 
                                                 << " Int Sin(Phi): "  << cosLUT[globalPhi - phi_quadrants[1]] << "\n";
          std::cout << "Int Phi: " << (float)globalPhi/binmap["globphi"] 
                                                 << " Int Cos(Phi): -" << (float)cosLUT[phi_quadrants[2]- 1 - globalPhi]/binmap["globphi"] 
                                                 << " Int Sin(Phi): "  << (float)cosLUT[globalPhi - phi_quadrants[1]]/binmap["globphi"] << "\n"; 
        }  
      }
      else if (globalPhi >= phi_quadrants[2] && globalPhi < phi_quadrants[3]){
        sumPx[Sector] = sumPx[Sector] - (pt*cosLUT[globalPhi - phi_quadrants[2]])/ binmap["globphi"]; 
        sumPy[Sector] = sumPy[Sector] - (pt*cosLUT[phi_quadrants[3] - 1 - globalPhi] )/ binmap["globphi"]; 

        if (debug == 2){
          std::cout << "Sector: " << Sector << " Quadrant: " << 3 << "\n";
          std::cout << "Int Phi: " << globalPhi 
                                                 << " Int Cos(Phi): -" << cosLUT[globalPhi - phi_quadrants[2]] 
                                                 << " Int Sin(Phi): -" << cosLUT[phi_quadrants[3] - 1- globalPhi] << "\n";
          std::cout << "Int Phi: " << (float)globalPhi/binmap["globphi"] 
                                                 << " Int Cos(Phi): -" << (float)cosLUT[globalPhi - phi_quadrants[2]]/binmap["globphi"] 
                                                 << " Int Sin(Phi): -" << (float)cosLUT[phi_quadrants[3] - 1 - globalPhi]/binmap["globphi"] << "\n";  
        }  
      }
      else if (globalPhi >= phi_quadrants[3] && globalPhi < phi_quadrants[4]){
        sumPx[Sector] = sumPx[Sector] + (pt*cosLUT[phi_quadrants[4] - 1 -  globalPhi] )/ binmap["globphi"];
        sumPy[Sector] = sumPy[Sector] - (pt*cosLUT[globalPhi - phi_quadrants[3]] )/ binmap["globphi"];

        if (debug == 2){
          std::cout << "Sector: " << Sector << " Quadrant: " << 4 << "\n";
          std::cout << "Int Phi: " << globalPhi 
                                                 << " Int Cos(Phi): "  << cosLUT[phi_quadrants[4] - 1 - globalPhi] 
                                                 << " Int Sin(Phi): -" << cosLUT[globalPhi -phi_quadrants[3]] << "\n";
          std::cout << "Int Phi: " << (float)globalPhi/binmap["globphi"] 
                                                 << " Int Cos(Phi): "  << (float)cosLUT[phi_quadrants[4] - 1 - globalPhi]/binmap["globphi"] 
                                                 << " Int Sin(Phi): -" << (float)cosLUT[globalPhi - phi_quadrants[3]]/binmap["globphi"] << "\n"; 
        }  
      }
    }
 
  } // end loop over tracks

  iEt GlobalPx = 0;
  iEt GlobalPy = 0;

  //Global Et sum with global phi bit shift to bring sums back to correct magnitude
  for (unsigned int i=0; i< N_sectors_; i++){
    GlobalPx = GlobalPx + sumPx[i];
    GlobalPy = GlobalPy + sumPy[i];
  }

  //Perform cordic sqrt, take x,y and converts to polar coordinate r/phi where r=sqrt(x**2+y**2) and phi = atan(y/x)
  EtMiss EtMiss = cordic_sqrt.to_polar(GlobalPx,GlobalPy);

  // Recentre phi
  int tempPhi = 0;

  if ((GlobalPx < 0) && (GlobalPy < 0))
    tempPhi =  EtMiss.Phi - phiScale_/2;
  else if ((GlobalPx >= 0) && (GlobalPy >= 0))
    tempPhi =  (EtMiss.Phi) + phiScale_/2;
  else if ((GlobalPx >= 0) && (GlobalPy < 0))
    tempPhi = EtMiss.Phi - phiScale_/2;
  else if ((GlobalPx < 0) && (GlobalPy >= 0))
    tempPhi = -EtMiss.Phi + 3*phiScale_/2;
  
  
  float EtScale = maxvalues["pt"]/binmap["pt"];
  float EtphiScale = (2*M_PI)/phiScale_;

  float converted_phi = (float)tempPhi*EtphiScale;

  if (debug == 4){
    std::cout << "====Sector Pt====" << "\n";
    for (auto i: sumPx){ std::cout << i*EtScale << "|";}
    std::cout << "\n";
    for (auto i: sumPy){ std::cout << i*EtScale << "|";}
    std::cout << "\n";
  
    std::cout << "====Global Pt====" << "\n";
    std::cout << GlobalPx << "|" << GlobalPy << "\n";
    std::cout << GlobalPx*EtScale << "|" << GlobalPy*EtScale << "\n";
  }

  if (debug == 6){
    std::cout << "====MET===" << "\n";
    std::cout << EtMiss.Et << "|" << EtMiss.Phi << "\n";
    std::cout << (EtMiss.Et)*EtScale << "|" << converted_phi << "\n";

    std::cout << "numTracks: " << num_tracks << "\n";
    std::cout << "quality Tracks: " << num_quality_tracks << "\n";
    std::cout << "assoc_tracks: " << num_assoc_tracks << "\n";
    std::cout << "========================================================" << "\n";
  }

  math::XYZTLorentzVector missingEt( -GlobalPx*EtScale, -GlobalPy*EtScale, 0, EtMiss.Et*EtScale);
  int ibx = 0;
  TkEtMiss L1EtObject( missingEt,
                       TkEtMiss::hwMET,
                       EtMiss.Et*EtScale,
                       converted_phi,
                       (unsigned int)num_assoc_tracks,
                       ibx);

  METCollection->push_back( L1EtObject );


  iEvent.put( std::move(METCollection), L1MetCollectionName);
} // end producer

void L1TrackerEtMissEmulatorProducer::beginJob() { }

void L1TrackerEtMissEmulatorProducer::endJob() { }

DEFINE_FWK_MODULE(L1TrackerEtMissEmulatorProducer);
