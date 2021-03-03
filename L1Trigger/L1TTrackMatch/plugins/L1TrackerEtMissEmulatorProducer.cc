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
  iEta      maxEta_; // in HWU
  ichi2XY   chi2rphidofMax_; // in HWU
  ichi2XY   chi2rzdofMax_; // in HWU
  ichi2bend bendChi2Max_; // in HWU
  iPt       minPt_; // in HWU
  iNstub    nStubsmin_;

  iglobPhi cosLUT[cosLUT_bins] = { 0 };// Cos LUT array
  iEta EtaRegionsLUT[N_etaregions+1]= { 0 }; //Various precomputed LUTs
  iZ0  DeltaZLUT[N_etaregions]= { 0 };
  iglobPhi phi_quadrants[N_quadrants]= { 0 };
  iglobPhi phi_shift[N_sectors+1]= { 0 };   

  int CordicSteps;
  int debug;
  bool cordicdebug = false;
  bool writeLUTs   = false;


  const edm::EDGetTokenT< TkPrimaryVertexCollection > pvToken_;
  const edm::EDGetTokenT<std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > > trackToken_;
};

//constructor//
L1TrackerEtMissEmulatorProducer::L1TrackerEtMissEmulatorProducer(const edm::ParameterSet& iConfig) :
pvToken_(consumes<TkPrimaryVertexCollection>(iConfig.getParameter<edm::InputTag>("L1VertexInputTag"))),
trackToken_(consumes< std::vector<TTTrack< Ref_Phase2TrackerDigi_> > > (iConfig.getParameter<edm::InputTag>("L1TrackInputTag")))
{
  // Input parameter cuts and convert to correct integer representations
  maxZ0_          = digitize_Signed <iZ0> ((float)iConfig.getParameter<double>("maxZ0"),-maxZ0, maxZ0, N_z0Bins);
  maxEta_         = digitize_Signed <iEta> ((float)iConfig.getParameter<double>("maxEta"), -maxEta, maxEta, N_etaBins);
  chi2rphidofMax_ = Chi_to_FWChi <ichi2XY,N_chi2XYbins> ((float)iConfig.getParameter<double>("chi2rphidofMax"),chi2Values);
  chi2rzdofMax_   = Chi_to_FWChi <ichi2XY,N_chi2XYbins> ((float)iConfig.getParameter<double>("chi2rzdofMax"),chi2Values);
  bendChi2Max_    = Chi_to_FWChi <ichi2bend,N_chi2bendbins> ((float)iConfig.getParameter<double>("bendChi2Max"),chi2bendValues);
  minPt_          = digitize_Signed <iPt> ((float)iConfig.getParameter<double>("minPt"),0,maxPt,N_ptBins);
  nStubsmin_      = (iNstub)iConfig.getParameter<int>("nStubsmin");

  CordicSteps     = (int)iConfig.getParameter<int>("nCordicSteps");
  debug           = (int)iConfig.getParameter<int>("Debug");
  writeLUTs       = (bool)iConfig.getParameter<bool>("WriteLUTs");


  if (debug == 5){cordicdebug=true;}


  // Compute LUTs
  FillCosLUT(cosLUT,max_LUT_phi);
  generate_phi_slice_LUTs <N_quadrants> (phi_quadrants);
  generate_phi_slice_LUTs <N_sectors+1> (phi_shift);
  generate_EtaRegions(EtaRegions,EtaRegionsLUT);
  generate_DeltaZBins(DeltaZ,DeltaZLUT);

  //Write Out LUTs
  if (writeLUTs){
    writeLUTout<iglobPhi,cosLUT_bins>    (cosLUT,"cos",",");
    writeLUTout<iglobPhi,N_quadrants>    (phi_quadrants,"phiquadrants",",");
    writeLUTout<iglobPhi,N_phiShift>     (phi_shift,"phishift",",");
    writeLUTout<iEta,N_etaregions+1>   (EtaRegionsLUT,"etaregions",",");
    writeLUTout<iZ0,N_etaregions> (DeltaZLUT,"dzbins",",");
    }


  if (debug == 1){
    edm::LogVerbatim("L1TkEtMissEmulator") << "=================================================\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << "Phi Quadrants: ";
    for (auto i: phi_quadrants){ edm::LogVerbatim("L1TkEtMissEmulator") << i << "|";}
    edm::LogVerbatim("L1TkEtMissEmulator") << "\n";

    edm::LogVerbatim("L1TkEtMissEmulator") << "Phi Sector Shift: ";
    for (auto i: phi_shift){ edm::LogVerbatim("L1TkEtMissEmulator") << i << "|";}
    edm::LogVerbatim("L1TkEtMissEmulator") << "\n";
    }

  
  produces<TkEtMissCollection>("L1TrackerEmuEtMiss");
}

L1TrackerEtMissEmulatorProducer::~L1TrackerEtMissEmulatorProducer() { }

void L1TrackerEtMissEmulatorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::unique_ptr<TkEtMissCollection> METCollection(new TkEtMissCollection);

  edm::Handle<TkPrimaryVertexCollection> L1VertexHandle;
  iEvent.getByToken(pvToken_,L1VertexHandle);

  edm::Handle<L1TTTrackCollectionType> L1TTTrackHandle;
  iEvent.getByToken(trackToken_, L1TTTrackHandle);

  Cordic cordic_sqrt(phiScale,N_ptBits, CordicSteps,cordicdebug,writeLUTs);


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
  iEt sumPx[N_sectors] = {  0,0,0,0,0,0,0,0,0 };
  iEt sumPy[N_sectors] = {  0,0,0,0,0,0,0,0,0 };

  // Track counters
  int num_tracks = 0;
  int num_assoc_tracks = 0;
  int num_quality_tracks = 0;

  iZ0 zVTX = digitize_Signed <iZ0> (L1VertexHandle->begin()->zvertex(), -maxZ0, maxZ0, N_z0Bins); //Convert vertex TODO use emulator

  for (const auto& track : *L1TTTrackHandle) {
    num_tracks++;
    L1TTTrackType& track_ref = const_cast< L1TTTrackType&>(track); //Non Const member functions in TTTrack_TrackWord 

    //Read in and convert track parameters to integer representation
    iPt pt   = digitize_Signed <iPt> (track_ref.momentum().perp(),0,maxPt,N_ptBins);
    iEta eta = digitize_Signed <iEta> (abs(track_ref.momentum().eta()),-maxEta, maxEta, N_etaBins);

    ichi2XY chi2rphidof = Chi_to_FWChi <ichi2XY,N_chi2XYbins> (track_ref.chi2XY(),chi2Values);   
    ichi2XY chi2rzdof   = Chi_to_FWChi <ichi2XY,N_chi2XYbins> (track_ref.chi2Z(),chi2Values); //TTTrack_TrackWord has no chi2Z getter function 
    ichi2bend bendChi2  = Chi_to_FWChi <ichi2bend,N_chi2bendbins> (track_ref.stubPtConsistency(),chi2bendValues);
    iNstub nstubs       = CountNStub(track_ref.unpack_hitPattern());

    unsigned int Sector  = track_ref.phiSector();
    // convert to local phi, functionality not in TTTrack_TrackWord word
    iPhi localPhi        = FloatPhi_to_iPhi(track_ref.phi(),Sector);
    // Convert to global phi
    iglobPhi globalPhi   = local_to_global(localPhi,phi_shift[Sector],phi_quadrants);

    iZ0  z0             = digitize_Signed <iZ0> (track_ref.z0(), -maxZ0, maxZ0, N_z0Bins);
    iZ0  deltaZ0        = 0;


    // Parameter cuts
    if (pt <= minPt_) continue;  
    if (z0 > maxZ0_) continue;
    if (eta > maxEta_) continue;
    if (chi2rphidof >= chi2rphidofMax_) continue;
    if (chi2rzdof >= chi2rzdofMax_) continue;
    if (bendChi2 >= bendChi2Max_) continue;
    if (nstubs < nStubsmin_) continue;

    num_quality_tracks++;
    // Temporary deltaZ to facilitate use of unsigned int for iZ0
    int temp = z0 - zVTX;
    iZ0 z_diff = abs(temp);
     
    // Track to vertex association
    if      ( eta>=EtaRegionsLUT[0] &&  eta< EtaRegionsLUT[1])  deltaZ0 = DeltaZLUT[0];
    else if ( eta>=EtaRegionsLUT[1] &&  eta< EtaRegionsLUT[2])  deltaZ0 = DeltaZLUT[1];
    else if ( eta>=EtaRegionsLUT[2] &&  eta< EtaRegionsLUT[3])  deltaZ0 = DeltaZLUT[2];
    else if ( eta>=EtaRegionsLUT[3] &&  eta< EtaRegionsLUT[4])  deltaZ0 = DeltaZLUT[3];
    else if ( eta>=EtaRegionsLUT[4] &&  eta< EtaRegionsLUT[5])  deltaZ0 = DeltaZLUT[4];
    else if ( eta>=EtaRegionsLUT[5] &&  eta<=EtaRegionsLUT[6])  deltaZ0 = DeltaZLUT[5];


    edm::LogVerbatim("L1TkEtMissEmulator") << "========================Z0 Debug=================================\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << "z0: " << z0 << " z vertex: " << zVTX << '\n';
    edm::LogVerbatim("L1TkEtMissEmulator") << "Zdiff: " << z_diff << " dZ bin: " << deltaZ0 << '\n';
    edm::LogVerbatim("L1TkEtMissEmulator") << "eta: " << eta << '\n';

    if ( z_diff <= deltaZ0 ) {
      num_assoc_tracks++;

      if (debug == 2){
        edm::LogVerbatim("L1TkEtMissEmulator") << "========================Phi Debug=================================\n";
        edm::LogVerbatim("L1TkEtMissEmulator") << "Sector Phi: " << localPhi << " Global Phi: " << globalPhi << "\n";
        edm::LogVerbatim("L1TkEtMissEmulator") << "Int Phi: "    << globalPhi << " Float Phi: " << track_ref.phi() 
                                               << " Float Cos(Phi): " << cos(track_ref.phi()) 
                                               << " Float Sin(Phi): " << sin(track_ref.phi()) << "\n";
      }
      
      // Split tracks in phi quadrants and access cosLUT, backwards iteration through cosLUT gives sin
      // Sum sector Et -ve when cos or sin phi are -ve
      if (globalPhi >= phi_quadrants[0] && globalPhi < phi_quadrants[1]){
        sumPx[Sector] = sumPx[Sector] + (pt*cosLUT[globalPhi] ); 
        sumPy[Sector] = sumPy[Sector] + (pt*cosLUT[phi_quadrants[1] - 1 - globalPhi]);

        if (debug == 2){
          edm::LogVerbatim("L1TkEtMissEmulator") << "Sector: " << Sector << " Quadrant: " << 1 << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator") << "Int Phi: " << globalPhi 
                                                 << " Int Cos(Phi): " << cosLUT[globalPhi] 
                                                 << " Int Sin(Phi): " << cosLUT[phi_quadrants[1]- 1 - globalPhi] << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator") << "Int Phi: " << (float)globalPhi/N_globPhiBins 
                                                 << " Int Cos(Phi): " << (float)cosLUT[globalPhi]/N_globPhiBins 
                                                 << " Int Sin(Phi): " << (float)cosLUT[phi_quadrants[1] - 1- globalPhi]/N_globPhiBins << "\n";  
        }   
      }
      else if (globalPhi >= phi_quadrants[1] && globalPhi < phi_quadrants[2]){
        sumPx[Sector] = sumPx[Sector] - (pt*cosLUT[phi_quadrants[2] - 1 - globalPhi] );  
        sumPy[Sector] = sumPy[Sector] + (pt*cosLUT[globalPhi - phi_quadrants[1]] );  

        if (debug == 2){
          edm::LogVerbatim("L1TkEtMissEmulator") << "Sector: " << Sector << " Quadrant: " << 2 << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator") << "Int Phi: " << globalPhi 
                                                 << " Int Cos(Phi): -" << cosLUT[phi_quadrants[2] - globalPhi] 
                                                 << " Int Sin(Phi): "  << cosLUT[globalPhi - phi_quadrants[1]] << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator") << "Int Phi: " << (float)globalPhi/N_globPhiBins 
                                                 << " Int Cos(Phi): -" << (float)cosLUT[phi_quadrants[2]- 1 - globalPhi]/N_globPhiBins 
                                                 << " Int Sin(Phi): "  << (float)cosLUT[globalPhi - phi_quadrants[1]]/N_globPhiBins << "\n"; 
        }  
      }
      else if (globalPhi >= phi_quadrants[2] && globalPhi < phi_quadrants[3]){
        sumPx[Sector] = sumPx[Sector] - (pt*cosLUT[globalPhi - phi_quadrants[2]]); 
        sumPy[Sector] = sumPy[Sector] - (pt*cosLUT[phi_quadrants[3] - 1 - globalPhi] ); 

        if (debug == 2){
          edm::LogVerbatim("L1TkEtMissEmulator") << "Sector: " << Sector << " Quadrant: " << 3 << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator") << "Int Phi: " << globalPhi 
                                                 << " Int Cos(Phi): -" << cosLUT[globalPhi - phi_quadrants[2]] 
                                                 << " Int Sin(Phi): -" << cosLUT[phi_quadrants[3] - 1- globalPhi] << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator") << "Int Phi: " << (float)globalPhi/N_globPhiBins 
                                                 << " Int Cos(Phi): -" << (float)cosLUT[globalPhi - phi_quadrants[2]]/N_globPhiBins 
                                                 << " Int Sin(Phi): -" << (float)cosLUT[phi_quadrants[3] - 1 - globalPhi]/N_globPhiBins << "\n";  
        }  
      }
      else if (globalPhi >= phi_quadrants[3] && globalPhi < phi_quadrants[4]){
        sumPx[Sector] = sumPx[Sector] + (pt*cosLUT[phi_quadrants[4] - 1 -  globalPhi] );
        sumPy[Sector] = sumPy[Sector] - (pt*cosLUT[globalPhi - phi_quadrants[3]] );

        if (debug == 2){
          edm::LogVerbatim("L1TkEtMissEmulator") << "Sector: " << Sector << " Quadrant: " << 4 << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator") << "Int Phi: " << globalPhi 
                                                 << " Int Cos(Phi): "  << cosLUT[phi_quadrants[4] - 1 - globalPhi] 
                                                 << " Int Sin(Phi): -" << cosLUT[globalPhi -phi_quadrants[3]] << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator") << "Int Phi: " << (float)globalPhi/N_globPhiBins 
                                                 << " Int Cos(Phi): "  << (float)cosLUT[phi_quadrants[4] - 1 - globalPhi]/N_globPhiBins 
                                                 << " Int Sin(Phi): -" << (float)cosLUT[globalPhi - phi_quadrants[3]]/N_globPhiBins << "\n"; 
        }  
      }
    }
 
  } // end loop over tracks

  
  iEt GlobalPx = 0;
  iEt GlobalPy = 0;

  //Global Et sum with global phi bit shift to bring sums back to correct magnitude
  for (int i=0; i< N_sectors; i++){
    GlobalPx = GlobalPx + sumPx[i]/ N_globPhiBins;
    GlobalPy = GlobalPy + sumPy[i]/ N_globPhiBins;
  }

  //Perform cordic sqrt, take x,y and converts to polar coordinate r/phi where r=sqrt(x**2+y**2) and phi = atan(y/x)
  EtMiss EtMiss_ = cordic_sqrt.to_polar(GlobalPx,GlobalPy);

  // Recentre phi
  
  int tempPhi = 0;

  if ((GlobalPx < 0) && (GlobalPy < 0))
    tempPhi =  EtMiss_.Phi - phiScale/2;
  else if ((GlobalPx >= 0) && (GlobalPy >= 0))
    tempPhi =  (EtMiss_.Phi) + phiScale/2;
  else if ((GlobalPx >= 0) && (GlobalPy < 0))
    tempPhi = EtMiss_.Phi - phiScale/2;
  else if ((GlobalPx < 0) && (GlobalPy >= 0))
    tempPhi = -EtMiss_.Phi + 3*phiScale/2;
  
  
  float EtScale = maxPt/N_ptBins;
  float EtphiScale = (2*M_PI)/phiScale;

  float converted_phi = (float)tempPhi*EtphiScale;

  if (debug == 4){
    edm::LogVerbatim("L1TkEtMissEmulator") << "====Sector Pt====" << "\n";
    for (auto i: sumPx){ edm::LogVerbatim("L1TkEtMissEmulator") << i*EtScale << "|";}
    edm::LogVerbatim("L1TkEtMissEmulator") << "\n";
    for (auto i: sumPy){ edm::LogVerbatim("L1TkEtMissEmulator") << i*EtScale << "|";}
    edm::LogVerbatim("L1TkEtMissEmulator") << "\n";
  
    edm::LogVerbatim("L1TkEtMissEmulator") << "====Global Pt====" << "\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << GlobalPx << "|" << GlobalPy << "\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << GlobalPx*EtScale << "|" << GlobalPy*EtScale << "\n";
  }

  if (debug == 6){
    edm::LogVerbatim("L1TkEtMissEmulator") << "====MET===" << "\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << EtMiss_.Et << "|" << EtMiss_.Phi << "\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << (EtMiss_.Et)*EtScale << "|" << converted_phi << "\n";

    edm::LogVerbatim("L1TkEtMissEmulator") << "numTracks: " << num_tracks << "\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << "quality Tracks: " << num_quality_tracks << "\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << "assoc_tracks: " << num_assoc_tracks << "\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << "========================================================" << "\n";
  }

  math::XYZTLorentzVector missingEt( -GlobalPx*EtScale, -GlobalPy*EtScale, 0, EtMiss_.Et*EtScale);
  int ibx = 0;
  TkEtMiss L1EtObject( missingEt,
                       TkEtMiss::hwMET,
                       EtMiss_.Et,
                       tempPhi,
                       num_assoc_tracks,
                       EtScale,
                       EtphiScale ,
                       ibx);
  L1EtObject.ihwtofloat();
  METCollection->push_back( L1EtObject );


  iEvent.put( std::move(METCollection), "L1TrackerEmuEtMiss");
} // end producer

void L1TrackerEtMissEmulatorProducer::beginJob() { }

void L1TrackerEtMissEmulatorProducer::endJob() { }

DEFINE_FWK_MODULE(L1TrackerEtMissEmulatorProducer);
