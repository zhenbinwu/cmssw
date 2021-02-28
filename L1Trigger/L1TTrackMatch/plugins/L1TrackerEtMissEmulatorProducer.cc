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
  ichi2XY   chi2rphidofMax_;
  ichi2XY   chi2rzdofMax_;
  ichi2bend bendChi2Max_;
  iPt       minPt_; // in HWU
  iNstub    nStubsmin_;

  iglobPhi cosLUT[LUT_bins];
  std::vector<iPhi> EtaRegionsLUT;
  std::vector<iZ0>  DeltaZLUT;
  std::vector<iglobPhi> phi_quadrants;
  std::vector<iglobPhi> phi_shift; 

  const edm::EDGetTokenT< TkPrimaryVertexCollection > pvToken_;
  const edm::EDGetTokenT<std::vector<TTTrack< Ref_Phase2TrackerDigi_ > > > trackToken_;
};

//constructor//
L1TrackerEtMissEmulatorProducer::L1TrackerEtMissEmulatorProducer(const edm::ParameterSet& iConfig) :
pvToken_(consumes<TkPrimaryVertexCollection>(iConfig.getParameter<edm::InputTag>("L1VertexInputTag"))),
trackToken_(consumes< std::vector<TTTrack< Ref_Phase2TrackerDigi_> > > (iConfig.getParameter<edm::InputTag>("L1TrackInputTag")))
{

  maxZ0_          = digitize_Signed <iZ0> ((float)iConfig.getParameter<double>("maxZ0"),-maxZ0, maxZ0, N_z0Bins);
  maxEta_         = digitize_Signed <iEta> ((float)iConfig.getParameter<double>("maxEta"), -maxEta, maxEta, N_etaBins);
  chi2rphidofMax_ = Chi_to_FWChi <ichi2XY,N_chi2XYbins> ((float)iConfig.getParameter<double>("chi2rphidofMax"),chi2Values);
  chi2rzdofMax_   = Chi_to_FWChi <ichi2XY,N_chi2XYbins> ((float)iConfig.getParameter<double>("chi2rzdofMax"),chi2Values);
  bendChi2Max_    = Chi_to_FWChi <ichi2bend,N_chi2bendbins> ((float)iConfig.getParameter<double>("bendChi2Max"),chi2bendValues);
  minPt_          = digitize_Signed <iPt> ((float)iConfig.getParameter<double>("minPt"),0,maxPt,N_ptBins);
  nStubsmin_      = (iNstub)iConfig.getParameter<int>("nStubsmin");

  FillCosLUT(cosLUT,max_LUT_phi);
  phi_quadrants = generate_phi_slice_LUTs(5);
  phi_shift = generate_phi_slice_LUTs(10);
  EtaRegionsLUT = generate_EtaRegions(EtaRegions);
  DeltaZLUT = generate_DeltaZBins(DeltaZ);

  

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

  Cordic cordic_sqrt(phiScale,N_ptBits, cordic_steps);


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

  iEt sumPx[9] = {0,0,0,0,0,0,0,0,0};
  iEt sumPy[9] = {0,0,0,0,0,0,0,0,0};
  int num_tracks = 0;
  int num_assoc_tracks = 0;
  int num_quality_tracks = 0;

  iZ0 zVTX = digitize_Signed <iZ0> (L1VertexHandle->begin()->zvertex(), -maxZ0, maxZ0, N_z0Bins);

  for (const auto& track : *L1TTTrackHandle) {
    num_tracks++;
    L1TTTrackType& track_ref = const_cast< L1TTTrackType&>(track); //Non Const member functions in TTTrack_TrackWord 

    iPt pt   = digitize_Signed <iPt> (track_ref.momentum().perp(),0,maxPt,N_ptBins);
    iEta eta = digitize_Signed <iEta> (abs(track_ref.momentum().eta()),-maxEta, maxEta, N_etaBins);

    ichi2XY chi2rphidof = Chi_to_FWChi <ichi2XY,N_chi2XYbins> (track_ref.chi2XY(),chi2Values);   
    ichi2XY chi2rzdof   = Chi_to_FWChi <ichi2XY,N_chi2XYbins> (track_ref.chi2Z(),chi2Values); //TTTrack_TrackWord has no chi2Z getter function 
    ichi2bend bendChi2  = Chi_to_FWChi <ichi2bend,N_chi2bendbins> (track_ref.stubPtConsistency(),chi2bendValues);
    iNstub nstubs       = CountNStub(track_ref.unpack_hitPattern());

    unsigned int Sector  = track_ref.phiSector();
    iPhi localPhi        = FloatPhi_to_iPhi(track_ref.phi(),Sector);
    iglobPhi globalPhi   = local_to_global(localPhi,phi_shift[Sector],phi_quadrants);

    iZ0  z0             = digitize_Signed <iZ0> (track_ref.z0(), -maxZ0, maxZ0, N_z0Bins);
    iZ0  deltaZ0        = 0;



    //std::cout << "====Phi Params====" << std::endl;
    //std::cout << Sector  << "|" << localPhi << "|" << globalPhi << "|" << track_ref.phi() << std::endl;
    
    //std::cout << "====Quality Params====" << std::endl;
    //std::cout << z0 << "|" << pt << "|" << eta << "|" << std::endl;
    //std::cout << track_ref.z0() << "|" << track_ref.momentum().perp() << "|" << track_ref.momentum().eta() << "|" << std::endl;
    
    
    
    if (pt <= minPt_) continue;  
    if (z0 > maxZ0_) continue;
    if (eta > maxEta_) continue;
    if (chi2rphidof >= chi2rphidofMax_) continue;
    if (chi2rzdof >= chi2rzdofMax_) continue;
    if (bendChi2 >= bendChi2Max_) continue;
    if (nstubs < nStubsmin_) continue;

    num_quality_tracks++;
    int temp = z0 - zVTX;
    iZ0 z_diff = abs(temp);
     

    if      ( eta>=EtaRegionsLUT[0] &&  eta< EtaRegionsLUT[1])  deltaZ0 = DeltaZLUT[0];
    else if ( eta>=EtaRegionsLUT[1] &&  eta< EtaRegionsLUT[2])  deltaZ0 = DeltaZLUT[1];
    else if ( eta>=EtaRegionsLUT[2] &&  eta< EtaRegionsLUT[3])  deltaZ0 = DeltaZLUT[2];
    else if ( eta>=EtaRegionsLUT[3] &&  eta< EtaRegionsLUT[4])  deltaZ0 = DeltaZLUT[3];
    else if ( eta>=EtaRegionsLUT[4] &&  eta< EtaRegionsLUT[5])  deltaZ0 = DeltaZLUT[4];
    else if ( eta>=EtaRegionsLUT[5] &&  eta<=EtaRegionsLUT[6])  deltaZ0 = DeltaZLUT[5];


    //std::cout << "Zdiff: " << z_diff << " dZ: " << deltaZ0 << std::endl;
    //std::cout << EtaRegions[0] << "|" << EtaRegions[1] << "|" << EtaRegions[2] << "|" << EtaRegions[3] << "|"
    //          << EtaRegions[4] << "|" << EtaRegions[5] << "|" << EtaRegions[6] << std::endl;
    //std::cout << DeltaZ[0] << "|" << DeltaZ[1] << "|" << DeltaZ[2] << "|" << DeltaZ[3] << "|"
    //          << DeltaZ[4] << "|" << DeltaZ[5] << "|" << std::endl;
    
    if ( z_diff <= deltaZ0) {
      num_assoc_tracks++;
      //std::cout << track_ref.momentum().perp() << "|" << track_ref.z0() << "|" << track_ref.momentum().eta() << "|" << track_ref.stubPtConsistency() << "|" << track_ref.chi2Z() << "|" << track_ref.chi2XY() << "|" << nstubs << std::endl;
      //std::cout << pt << "|" << bendChi2 << "|" << chi2rzdof << "|" << chi2rphidof << "|" << std::endl;
      //sumPx[Sector] = sumPx[Sector] + pt*cos(track_ref.phi());
      //sumPy[Sector] = sumPy[Sector] + pt*sin(track_ref.phi());
      //std::cout << "===============================================================" << std::endl;
      //std::cout << "Sector Phi: " << localPhi << " Global Phi: " << globalPhi << std::endl;
      //std::cout << "Int Phi: " << globalPhi << " Float Phi: " << track_ref.phi() << " Float Cos(Phi): " << cos(track_ref.phi()) << " Float Sin(Phi): " << sin(track_ref.phi()) << std::endl;
      //std::cout << phi_quadrants[0] <<  "|" << phi_quadrants[1] << "|" << phi_quadrants[2] << "|" << phi_quadrants[3] << "|" << phi_quadrants[4] << std::endl;
      
      
      if (globalPhi >= phi_quadrants[0] && globalPhi < phi_quadrants[1]){
        sumPx[Sector] = sumPx[Sector] + (pt*cosLUT[globalPhi] ) ;  
        sumPy[Sector] = sumPy[Sector] + (pt*cosLUT[phi_quadrants[1] - 1 - globalPhi]) ;
        //std::cout << "Int Phi: " << globalPhi << " Int Cos(Phi): " << cosLUT[globalPhi] << " Int Sin(Phi): " << cosLUT[phi_quadrants[1] - 1 - globalPhi] << std::endl;
        //std::cout << "Int Phi: " << (float)globalPhi/N_globPhiBins << " Int Cos(Phi): " << (float)cosLUT[globalPhi]/N_globPhiBins << " Int Sin(Phi): " << (float)cosLUT[phi_quadrants[1] - 1 - globalPhi]/N_globPhiBins << std::endl;     
      }
      else if (globalPhi >= phi_quadrants[1] && globalPhi < phi_quadrants[2]){
        sumPx[Sector] = sumPx[Sector] - (pt*cosLUT[phi_quadrants[2] - 1 - globalPhi] );  
        sumPy[Sector] = sumPy[Sector] + (pt*cosLUT[globalPhi - phi_quadrants[1]] );  
        //std::cout << "Int Phi: " << globalPhi << " Int Cos(Phi): -" << cosLUT[phi_quadrants[2] - 1 - globalPhi] << " Int Sin(Phi): " << cosLUT[globalPhi - phi_quadrants[1]] << std::endl;
        //std::cout << "Int Phi: " << (float)globalPhi/N_globPhiBins << " Int Cos(Phi): -" << (float)cosLUT[phi_quadrants[2] - 1 - globalPhi]/N_globPhiBins << " Int Sin(Phi): " << (float)cosLUT[globalPhi - phi_quadrants[1]]/N_globPhiBins << std::endl;     
      }
      else if (globalPhi >= phi_quadrants[2] && globalPhi < phi_quadrants[3]){
        sumPx[Sector] = sumPx[Sector] - (pt*cosLUT[globalPhi - phi_quadrants[2]]);  
        sumPy[Sector] = sumPy[Sector] - (pt*cosLUT[phi_quadrants[3] - 1 - globalPhi] ); 
        //std::cout << "Int Phi: " << globalPhi << " Int Cos(Phi): -" << cosLUT[globalPhi - phi_quadrants[2]] << " Int Sin(Phi): " << -cosLUT[phi_quadrants[3] - 1 - globalPhi] << std::endl;
        //std::cout << "Int Phi: " << (float)globalPhi/N_globPhiBins << " Int Cos(Phi): -" << (float)cosLUT[globalPhi - phi_quadrants[2]]/N_globPhiBins << " Int Sin(Phi): " << -(float)cosLUT[phi_quadrants[3] - 1 - globalPhi]/N_globPhiBins << std::endl;     
      }
      else if (globalPhi >= phi_quadrants[3] && globalPhi < phi_quadrants[4]){
        sumPx[Sector] = sumPx[Sector] + (pt*cosLUT[phi_quadrants[4] - 1 - globalPhi] );  
        sumPy[Sector] = sumPy[Sector] - (pt*cosLUT[globalPhi - phi_quadrants[3]] );  
        //std::cout << "Int Phi: " << globalPhi << " Int Cos(Phi): " << cosLUT[phi_quadrants[4] - 1 - globalPhi] << " Int Sin(Phi): -" << cosLUT[globalPhi - phi_quadrants[3]] << std::endl;
        //std::cout << "Int Phi: " << (float)globalPhi/N_globPhiBins << " Int Cos(Phi): " << (float)cosLUT[phi_quadrants[4] - 1 - globalPhi]/N_globPhiBins << " Int Sin(Phi): -" << (float)cosLUT[globalPhi - phi_quadrants[3]]/N_globPhiBins << std::endl;     
      }
    }
 
  } // end loop over tracks
   
  iEt GlobalPx = 0;
  iEt GlobalPy = 0;


  for (int i=0; i< 9; i++){
    GlobalPx = GlobalPx + sumPx[i]/ pow(2,N_globphiBits);
    GlobalPy = GlobalPy + sumPy[i]/ pow(2,N_globphiBits);
  }

  EtMiss EtMiss_ = cordic_sqrt.to_polar(GlobalPx,GlobalPy);

  //Convert from Emulator quantities back to float to compare to simulation
  
  float converted_phi = EtMiss_.Phi*(2*M_PI)/phiScale;

  if ((EtMiss_.Phi > phiScale/2) && (GlobalPx >= 0) && (GlobalPy < 0))
    converted_phi =  (converted_phi - 2*M_PI);
  else if ((EtMiss_.Phi > phiScale/2) && (GlobalPx <  0) && (GlobalPy >= 0))
    converted_phi =  -(converted_phi - 2*M_PI);
  else if (EtMiss_.Phi > phiScale/2)
    converted_phi = (converted_phi - 2*M_PI);
  
  float EtScale = maxPt/N_ptBins;

  /*
  std::cout << "====Sector Pt====" << std::endl;

  std::cout << (sumPx[0] *EtScale) << "|" << (sumPx[1] *EtScale) << "|" 
            << (sumPx[2] *EtScale) << "|" << (sumPx[3] *EtScale) << "|" 
            << (sumPx[4] *EtScale) << "|" << (sumPx[5] *EtScale) << "|" 
            << (sumPx[6] *EtScale) << "|" << (sumPx[7] *EtScale) << "|" 
            << (sumPx[8] *EtScale) << std::endl;
  std::cout << (sumPy[0] *EtScale) << "|" << (sumPy[1] *EtScale) << "|" 
            << (sumPy[2] *EtScale) << "|" << (sumPy[3] *EtScale) << "|" 
            << (sumPy[4] *EtScale) << "|" << (sumPy[5] *EtScale) << "|" 
            << (sumPy[6] *EtScale) << "|" << (sumPy[7] *EtScale) << "|" 
            << (sumPy[8] *EtScale) << std::endl;
  */
  std::cout << "====Global Pt====" << std::endl;
  std::cout << GlobalPx << "|" << GlobalPy << std::endl;
  std::cout << GlobalPx*EtScale << "|" << GlobalPy*EtScale << std::endl;

  std::cout << "====MET===" << std::endl;
  std::cout << EtMiss_.Et << "|" << EtMiss_.Phi << std::endl;
  std::cout << (EtMiss_.Et)*EtScale << "|" << converted_phi << std::endl;

  std::cout << "numTracks: " << num_tracks << std::endl;
  std::cout << "quality Tracks: " << num_quality_tracks << std::endl;
  std::cout << "assoc_tracks: " << num_assoc_tracks << std::endl;
  std::cout << "========================================================" << std::endl;
  


  

  math::XYZTLorentzVector missingEt( -GlobalPx*EtScale, -GlobalPy*EtScale, 0, EtMiss_.Et*EtScale);

  METCollection->push_back( TkEtMiss( missingEt,
    TkEtMiss::kMET,
    converted_phi,
    0,
    0,
    0 )
  );


  iEvent.put( std::move(METCollection), "L1TrackerEmuEtMiss");
} // end producer

void L1TrackerEtMissEmulatorProducer::beginJob() { }

void L1TrackerEtMissEmulatorProducer::endJob() { }

DEFINE_FWK_MODULE(L1TrackerEtMissEmulatorProducer);
