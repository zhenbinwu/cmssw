// Original Author:  Christopher Brown
//         Created:  19 Feb 2021

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
#include "DataFormats/L1Trigger/interface/EtSum.h"
#include "DataFormats/L1TCorrelator/interface/TkPrimaryVertex.h"

#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuAlgo.h"
#include "L1Trigger/L1TTrackMatch/interface/Cordic.h"
#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuTrackTransform.h"

#include <numeric>

using namespace L1TkEtMissEmuAlgo;
using namespace l1t;

class L1TrackerEtMissEmulatorProducer : public edm::stream::EDProducer<> {
public:
  typedef TTTrack<Ref_Phase2TrackerDigi_> L1TTTrackType;
  typedef std::vector<L1TTTrackType> L1TTTrackCollectionType;

  explicit L1TrackerEtMissEmulatorProducer(const edm::ParameterSet&);
  ~L1TrackerEtMissEmulatorProducer() override;

private:
  virtual void beginJob();
  void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob();

  // ----------member data ---------------------------
  InternalEtWord EtCuts;

  std::vector<iglobPhi> cosLUT;     // Cos LUT array
  std::vector<iEta> EtaRegionsLUT;  //Various precomputed LUTs
  std::vector<iZ0> DeltaZLUT;
  std::vector<iglobPhi> phi_quadrants;
  std::vector<iglobPhi> phi_shifts;

  int CordicSteps;
  int debug;
  bool cordicdebug = false;
  bool writeLUTs = false;

  L1TkEtMissEmuTrackTransform TrackTransform;

  std::string L1MetCollectionName;

  const edm::EDGetTokenT<TkPrimaryVertexCollection> pvToken_;
  const edm::EDGetTokenT<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>> trackToken_;
};

//constructor//
L1TrackerEtMissEmulatorProducer::L1TrackerEtMissEmulatorProducer(const edm::ParameterSet& iConfig)
    : pvToken_(consumes<TkPrimaryVertexCollection>(iConfig.getParameter<edm::InputTag>("L1VertexInputTag"))),
      trackToken_(consumes<std::vector<TTTrack<Ref_Phase2TrackerDigi_>>>(
          iConfig.getParameter<edm::InputTag>("L1TrackInputTag"))) {
  // Input parameter cuts and convert to correct integer representations
  TrackTransform.GenerateLUTs();
  EtCuts = TrackTransform.TransfromCuts(iConfig);
  phi_quadrants = TrackTransform.get_phi_quad();
  phi_shifts = TrackTransform.get_phi_shift();

  CordicSteps = (int)iConfig.getParameter<int>("nCordicSteps");
  debug = (int)iConfig.getParameter<int>("Debug");
  writeLUTs = (bool)iConfig.getParameter<bool>("WriteLUTs");

  L1MetCollectionName = (std::string)iConfig.getParameter<std::string>("L1MetCollectionName");

  if (debug == 5) {
    cordicdebug = true;
  }

  //To have same bin spacing between 0 and pi/2 as between original phi granularity
  int cosLUT_bins = ceil((L1TkEtMissEmuAlgo::max_LUT_phi_ * (L1TkEtMissEmuAlgo::N_globPhiBins_ - 1)) /
                         (2 * L1TkEtMissEmuAlgo::max_TWord_Phi_));

  // Compute LUTs
  cosLUT = FillCosLUT(cosLUT_bins);
  EtaRegionsLUT = generate_EtaRegions();
  DeltaZLUT = generate_DeltaZBins();

  //Write Out LUTs
  if (writeLUTs) {
    writevectorLUTout<iglobPhi>(cosLUT, "cos", ",");
    writevectorLUTout<iglobPhi>(phi_quadrants, "phiquadrants", ",");
    writevectorLUTout<iglobPhi>(phi_shifts, "phishift", ",");
    writevectorLUTout<iEta>(EtaRegionsLUT, "etaregions", ",");
    writevectorLUTout<iZ0>(DeltaZLUT, "dzbins", ",");
    std::vector<unsigned int> cutlut = {(unsigned int)EtCuts.pt,
                                        (unsigned int)EtCuts.z0,
                                        (unsigned int)EtCuts.eta,
                                        (unsigned int)EtCuts.chi2rphidof,
                                        (unsigned int)EtCuts.chi2rzdof,
                                        (unsigned int)EtCuts.bendChi2,
                                        (unsigned int)EtCuts.nstubs};
    writevectorLUTout<unsigned int>(cutlut, "cuts", ",");
  }

  produces<std::vector<EtSum>>(L1MetCollectionName);
}

L1TrackerEtMissEmulatorProducer::~L1TrackerEtMissEmulatorProducer() {}

void L1TrackerEtMissEmulatorProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  std::unique_ptr<std::vector<l1t::EtSum>> METCollection(new std::vector<l1t::EtSum>(0));

  edm::Handle<TkPrimaryVertexCollection> L1VertexHandle;
  iEvent.getByToken(pvToken_, L1VertexHandle);

  edm::Handle<L1TTTrackCollectionType> L1TTTrackHandle;
  iEvent.getByToken(trackToken_, L1TTTrackHandle);

  //Initialize cordic class
  Cordic cordicSqrt(
      L1TkEtMissEmuAlgo::N_METphibins_, L1TkEtMissEmuAlgo::N_METbits_, CordicSteps, cordicdebug, writeLUTs);

  if (!L1VertexHandle.isValid()) {
    LogError("L1TrackerEtMissEmulatorProducer")
        << "\nWarning: TkPrimaryVertexCollection not found in the event. Exit\n";
    return;
  }

  if (!L1TTTrackHandle.isValid()) {
    LogError("L1TrackerEtMissEmulatorProducer") << "\nWarning: L1TTTrackCollection not found in the event. Exit\n";
    return;
  }

  // Initialize sector sums, need 0 initialization in case a sector has no tracks
  iEt sumPx[L1TkEtMissEmuAlgo::N_sectors_] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  iEt sumPy[L1TkEtMissEmuAlgo::N_sectors_] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

  // Track counters
  int num_tracks{0};
  int num_assoc_tracks{0};
  int num_quality_tracks{0};

  float zVTX = L1VertexHandle->begin()->zvertex();

  for (const auto& track : *L1TTTrackHandle) {
    num_tracks++;
    L1TTTrackType& track_ref = const_cast<L1TTTrackType&>(track);  //Non Const member functions in TTTrack_TrackWord

    InternalEtWord EtTrack = TrackTransform.TransformTrack(track_ref, zVTX);
    iZ0 deltaZ0 = 0;

    // Parameter cuts
    if (EtTrack.pt <= EtCuts.pt)
      continue;
    if (EtTrack.z0 > EtCuts.z0)
      continue;
    if (EtTrack.eta > EtCuts.eta)
      continue;
    if (EtTrack.chi2rphidof >= EtCuts.chi2rphidof)
      continue;
    if (EtTrack.chi2rzdof >= EtCuts.chi2rphidof)
      continue;
    if (EtTrack.bendChi2 >= EtCuts.bendChi2)
      continue;
    if (EtTrack.nstubs < EtCuts.nstubs)
      continue;

    num_quality_tracks++;
    // Temporary deltaZ to facilitate use of unsigned int for iZ0
    int temp = EtTrack.z0 - EtTrack.pV;
    iZ0 z_diff = abs(temp);

    // Track to vertex association, adaptive z window based on eta region
    for (unsigned int reg = 0; reg < L1TkEtMissEmuAlgo::N_etaregions_; reg++) {
      if (EtTrack.eta >= EtaRegionsLUT[reg] && EtTrack.eta < EtaRegionsLUT[reg + 1]) {
        deltaZ0 = DeltaZLUT[reg];
        break;
      }
    }

    if (z_diff <= deltaZ0) {
      num_assoc_tracks++;

      if (debug == 2) {
        edm::LogVerbatim("L1TkEtMissEmulator")
            << "========================Phi Debug=================================\n";
        edm::LogVerbatim("L1TkEtMissEmulator") << "Int pT: " << EtTrack.pt << "\n";
        edm::LogVerbatim("L1TkEtMissEmulator")
            << "Int Phi: " << EtTrack.globalPhi << " Float Phi: " << EtTrack.phi
            << " Actual Float Cos(Phi): " << cos(EtTrack.phi) << " Actual Float Sin(Phi): " << sin(EtTrack.phi) << "\n";
      }

      // Split tracks in phi quadrants and access cosLUT, backwards iteration through cosLUT gives sin
      // Sum sector Et -ve when cos or sin phi are -ve
      if (EtTrack.globalPhi >= phi_quadrants[0] && EtTrack.globalPhi < phi_quadrants[1]) {
        sumPx[EtTrack.Sector] =
            sumPx[EtTrack.Sector] + (EtTrack.pt * cosLUT[EtTrack.globalPhi]) / L1TkEtMissEmuAlgo::N_globPhiBins_;
        sumPy[EtTrack.Sector] =
            sumPy[EtTrack.Sector] +
            (EtTrack.pt * cosLUT[phi_quadrants[1] - 1 - EtTrack.globalPhi]) / L1TkEtMissEmuAlgo::N_globPhiBins_;

        if (debug == 2) {
          edm::LogVerbatim("L1TkEtMissEmulator") << "Sector: " << EtTrack.Sector << " Quadrant: " << 1 << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator")
              << "Int Phi: " << EtTrack.globalPhi << " Int Cos(Phi): " << cosLUT[EtTrack.globalPhi]
              << " Int Sin(Phi): " << cosLUT[phi_quadrants[1] - 1 - EtTrack.globalPhi] << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator")
              << "Float Phi: " << (float)EtTrack.globalPhi / L1TkEtMissEmuAlgo::N_globPhiBins_
              << " Float Cos(Phi): " << (float)cosLUT[EtTrack.globalPhi] / L1TkEtMissEmuAlgo::N_globPhiBins_
              << " Float Sin(Phi): "
              << (float)cosLUT[phi_quadrants[1] - 1 - EtTrack.globalPhi] / L1TkEtMissEmuAlgo::N_globPhiBins_ << "\n";
        }
      } else if (EtTrack.globalPhi >= phi_quadrants[1] && EtTrack.globalPhi < phi_quadrants[2]) {
        sumPx[EtTrack.Sector] =
            sumPx[EtTrack.Sector] -
            (EtTrack.pt * cosLUT[phi_quadrants[2] - 1 - EtTrack.globalPhi]) / L1TkEtMissEmuAlgo::N_globPhiBins_;
        sumPy[EtTrack.Sector] = sumPy[EtTrack.Sector] + (EtTrack.pt * cosLUT[EtTrack.globalPhi - phi_quadrants[1]]) /
                                                            L1TkEtMissEmuAlgo::N_globPhiBins_;

        if (debug == 2) {
          edm::LogVerbatim("L1TkEtMissEmulator") << "Sector: " << EtTrack.Sector << " Quadrant: " << 2 << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator")
              << "Int Phi: " << EtTrack.globalPhi << " Int Cos(Phi): -" << cosLUT[phi_quadrants[2] - EtTrack.globalPhi]
              << " Int Sin(Phi): " << cosLUT[EtTrack.globalPhi - phi_quadrants[1]] << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator")
              << "Float Phi: " << (float)EtTrack.globalPhi / L1TkEtMissEmuAlgo::N_globPhiBins_ << " Float Cos(Phi): -"
              << (float)cosLUT[phi_quadrants[2] - 1 - EtTrack.globalPhi] / L1TkEtMissEmuAlgo::N_globPhiBins_
              << " Float Sin(Phi): "
              << (float)cosLUT[EtTrack.globalPhi - phi_quadrants[1]] / L1TkEtMissEmuAlgo::N_globPhiBins_ << "\n";
        }
      } else if (EtTrack.globalPhi >= phi_quadrants[2] && EtTrack.globalPhi < phi_quadrants[3]) {
        sumPx[EtTrack.Sector] = sumPx[EtTrack.Sector] - (EtTrack.pt * cosLUT[EtTrack.globalPhi - phi_quadrants[2]]) /
                                                            L1TkEtMissEmuAlgo::N_globPhiBins_;
        sumPy[EtTrack.Sector] =
            sumPy[EtTrack.Sector] -
            (EtTrack.pt * cosLUT[phi_quadrants[3] - 1 - EtTrack.globalPhi]) / L1TkEtMissEmuAlgo::N_globPhiBins_;

        if (debug == 2) {
          edm::LogVerbatim("L1TkEtMissEmulator") << "Sector: " << EtTrack.Sector << " Quadrant: " << 3 << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator")
              << "Int Phi: " << EtTrack.globalPhi << " Int Cos(Phi): -" << cosLUT[EtTrack.globalPhi - phi_quadrants[2]]
              << " Int Sin(Phi): -" << cosLUT[phi_quadrants[3] - 1 - EtTrack.globalPhi] << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator")
              << "Float Phi: " << (float)EtTrack.globalPhi / L1TkEtMissEmuAlgo::N_globPhiBins_ << " Float Cos(Phi): -"
              << (float)cosLUT[EtTrack.globalPhi - phi_quadrants[2]] / L1TkEtMissEmuAlgo::N_globPhiBins_
              << " Float Sin(Phi): -"
              << (float)cosLUT[phi_quadrants[3] - 1 - EtTrack.globalPhi] / L1TkEtMissEmuAlgo::N_globPhiBins_ << "\n";
        }

      } else if (EtTrack.globalPhi >= phi_quadrants[3] && EtTrack.globalPhi < phi_quadrants[4]) {
        sumPx[EtTrack.Sector] =
            sumPx[EtTrack.Sector] +
            (EtTrack.pt * cosLUT[phi_quadrants[4] - 1 - EtTrack.globalPhi]) / L1TkEtMissEmuAlgo::N_globPhiBins_;
        sumPy[EtTrack.Sector] = sumPy[EtTrack.Sector] - (EtTrack.pt * cosLUT[EtTrack.globalPhi - phi_quadrants[3]]) /
                                                            L1TkEtMissEmuAlgo::N_globPhiBins_;

        if (debug == 2) {
          edm::LogVerbatim("L1TkEtMissEmulator") << "Sector: " << EtTrack.Sector << " Quadrant: " << 4 << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator")
              << "Int Phi: " << EtTrack.globalPhi
              << " Int Cos(Phi): " << cosLUT[phi_quadrants[4] - 1 - EtTrack.globalPhi] << " Int Sin(Phi): -"
              << cosLUT[EtTrack.globalPhi - phi_quadrants[3]] << "\n";
          edm::LogVerbatim("L1TkEtMissEmulator")
              << "Float Phi: " << (float)EtTrack.globalPhi / L1TkEtMissEmuAlgo::N_globPhiBins_ << " Float Cos(Phi): "
              << (float)cosLUT[phi_quadrants[4] - 1 - EtTrack.globalPhi] / L1TkEtMissEmuAlgo::N_globPhiBins_
              << " Float Sin(Phi): -"
              << (float)cosLUT[EtTrack.globalPhi - phi_quadrants[3]] / L1TkEtMissEmuAlgo::N_globPhiBins_ << "\n";
        }
      }
    }

  }  // end loop over tracks

  iEt GlobalPx = 0;
  iEt GlobalPy = 0;

  //Global Et sum
  for (unsigned int i = 0; i < L1TkEtMissEmuAlgo::N_sectors_; i++) {
    GlobalPx = GlobalPx + sumPx[i];
    GlobalPy = GlobalPy + sumPy[i];
  }

  //Perform cordic sqrt, take x,y and converts to polar coordinate r,phi where r=sqrt(x**2+y**2) and phi = atan(y/x)
  EtMiss EtMiss = cordicSqrt.toPolar(GlobalPx, GlobalPy);

  // Recentre phi
  iMETphi tempPhi = 0;

  if ((GlobalPx < 0) && (GlobalPy < 0))
    tempPhi = EtMiss.Phi - L1TkEtMissEmuAlgo::N_METphibins_ / 2;
  else if ((GlobalPx >= 0) && (GlobalPy >= 0))
    tempPhi = (EtMiss.Phi) + L1TkEtMissEmuAlgo::N_METphibins_ / 2;
  else if ((GlobalPx >= 0) && (GlobalPy < 0))
    tempPhi = EtMiss.Phi - L1TkEtMissEmuAlgo::N_METphibins_ / 2;
  else if ((GlobalPx < 0) && (GlobalPy >= 0))
    tempPhi = -EtMiss.Phi + 3 * L1TkEtMissEmuAlgo::N_METphibins_ / 2;

  if (debug == 6) {
    float EtScale = (L1TkEtMissEmuAlgo::max_METWord_ET_ / L1TkEtMissEmuAlgo::N_METbins_);
    float EtphiScale = (L1TkEtMissEmuAlgo::max_METWord_Phi_ / L1TkEtMissEmuAlgo::N_METphibins_);
    edm::LogVerbatim("L1TkEtMissEmulator") << "====Sector Pt====\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << "Px: ";
    for (auto i : sumPx) {
      edm::LogVerbatim("L1TkEtMissEmulator") << i * EtScale << "|";
    }
    edm::LogVerbatim("L1TkEtMissEmulator") << "\n Py: ";
    for (auto i : sumPy) {
      edm::LogVerbatim("L1TkEtMissEmulator") << i * EtScale << "|";
    }
    edm::LogVerbatim("L1TkEtMissEmulator") << "\n";

    edm::LogVerbatim("L1TkEtMissEmulator") << "====Global Pt====\n";
    edm::LogVerbatim("L1TkEtMissEmulator")
        << "Integer Global Px: " << GlobalPx << "| Integer Global Py: " << GlobalPy << "\n";
    edm::LogVerbatim("L1TkEtMissEmulator")
        << "Float Global Px: " << GlobalPx * EtScale << "| Float Global Py: " << GlobalPy * EtScale << "\n";
  }

  if (debug == 6) {
    float EtScale = (L1TkEtMissEmuAlgo::max_METWord_ET_ / L1TkEtMissEmuAlgo::N_METbins_);
    float EtphiScale = (L1TkEtMissEmuAlgo::max_METWord_Phi_ / L1TkEtMissEmuAlgo::N_METphibins_);
    edm::LogVerbatim("L1TkEtMissEmulator") << "====MET===\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << "Integer MET: " << EtMiss.Et << "| Integer MET phi" << EtMiss.Phi << "\n";
    edm::LogVerbatim("L1TkEtMissEmulator")
        << "Float MET: " << (EtMiss.Et) * EtScale << "| Float MET phi" << (float)tempPhi * EtphiScale - M_PI << "\n";

    edm::LogVerbatim("L1TkEtMissEmulator") << "# Intial Tracks: " << num_tracks << "\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << "# Tracks after Quality Cuts: " << num_quality_tracks << "\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << "# Tracks Associated to Vertex: " << num_assoc_tracks << "\n";
    edm::LogVerbatim("L1TkEtMissEmulator") << "========================================================\n";
  }

  math::XYZTLorentzVector missingEt(-GlobalPx, -GlobalPy, 0, EtMiss.Et);
  EtSum L1EtSum(missingEt, EtSum::EtSumType::kMissingEt, (int)EtMiss.Et, 0, (int)tempPhi, (int)num_assoc_tracks);

  METCollection->push_back(L1EtSum);

  iEvent.put(std::move(METCollection), L1MetCollectionName);
}  // end producer

void L1TrackerEtMissEmulatorProducer::beginJob() {}

void L1TrackerEtMissEmulatorProducer::endJob() {}

DEFINE_FWK_MODULE(L1TrackerEtMissEmulatorProducer);
