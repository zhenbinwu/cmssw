/******* \class DTRunConditionVarClient *******
 *
 * Description:
 *  
 *  detailed description
 *
 * \author : Paolo Bellan, Antonio Branca
 * $date   : 23/09/2011 15:42:04 CET $
 *
 * Modification:
 *
 *  threadsafe version (//-) oct/nov 2014 - WATWanAbdullah -ncpp-um-my
 *
 *
 *********************************/

#include "DQM/DTMonitorClient/src/DTRunConditionVarClient.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"

#include <cstdio>
#include <sstream>
#include <cmath>

using namespace edm;
using namespace std;

DTRunConditionVarClient::DTRunConditionVarClient(const ParameterSet& pSet) {
  LogVerbatim("DTDQM|DTMonitorClient|DTRunConditionVarClient") << "DTRunConditionVarClient: Constructor called";

  minRangeVDrift = pSet.getUntrackedParameter<double>("minRangeVDrift");
  maxRangeVDrift = pSet.getUntrackedParameter<double>("maxRangeVDrift");
  minRangeT0 = pSet.getUntrackedParameter<double>("minRangeT0");
  maxRangeT0 = pSet.getUntrackedParameter<double>("maxRangeT0");

  maxGoodVDriftDev = pSet.getUntrackedParameter<double>("maxGoodVDriftDev");
  minBadVDriftDev = pSet.getUntrackedParameter<double>("minBadVDriftDev");
  maxGoodT0 = pSet.getUntrackedParameter<double>("maxGoodT0");
  minBadT0 = pSet.getUntrackedParameter<double>("minBadT0");

  maxGoodVDriftSigma = pSet.getUntrackedParameter<double>("maxGoodVDriftSigma");
  minBadVDriftSigma = pSet.getUntrackedParameter<double>("minBadVDriftSigma");
  maxGoodT0Sigma = pSet.getUntrackedParameter<double>("maxGoodT0Sigma");
  minBadT0Sigma = pSet.getUntrackedParameter<double>("minBadT0Sigma");

  readLegacyVDriftDB = pSet.getParameter<bool>("readLegacyVDriftDB");

  nevents = 0;

  bookingdone = false;

  if (readLegacyVDriftDB) {
    mTimeMapToken_ = esConsumes<edm::Transition::BeginRun>();
  } else {
    vDriftToken_ = esConsumes<edm::Transition::BeginRun>();
  }
}

DTRunConditionVarClient::~DTRunConditionVarClient() {
  LogVerbatim("DTDQM|DTMonitorClient|DTRunConditionVarClient") << "DTRunConditionVarClient: Destructor called";
}

void DTRunConditionVarClient::beginRun(const Run& run, const EventSetup& context) {
  LogTrace("DTDQM|DTMonitorClient|DTResolutionAnalysisTest") << "[DTRunConditionVarClient]: BeginRun";
  // Get the map of vdrift from the setup
  if (readLegacyVDriftDB) {
    mTimeMap_ = &context.getData(mTimeMapToken_);
    vDriftMap_ = nullptr;
  } else {
    vDriftMap_ = &context.getData(vDriftToken_);
    mTimeMap_ = nullptr;
    // Consistency check: no parametrization is implemented for the time being
    int version = vDriftMap_->version();
    if (version != 1) {
      throw cms::Exception("Configuration") << "only version 1 is presently supported for VDriftDB";
    }
  }
}

void DTRunConditionVarClient::dqmEndLuminosityBlock(DQMStore::IBooker& ibooker,
                                                    DQMStore::IGetter& igetter,
                                                    edm::LuminosityBlock const& lumiSeg,
                                                    edm::EventSetup const& context) {}

void DTRunConditionVarClient::dqmEndJob(DQMStore::IBooker& ibooker, DQMStore::IGetter& igetter) {
  LogVerbatim("DTDQM|DTMonitorClient|DTRunConditionVarClient") << "DTRunConditionVarClient: end job";

  ibooker.setCurrentFolder("DT/02-Segments");

  glbVDriftSummary =
      ibooker.book2D("VDriftGlbSummary", "# of MBs with good mean and good sigma of vDrift", 12, 1, 13, 5, -2, 3);
  glbT0Summary = ibooker.book2D("T0GlbSummary", "# of MBs with good mean and good sigma of t0", 12, 1, 13, 5, -2, 3);

  ibooker.setCurrentFolder("DT/02-Segments/02-MeanVDrift");

  summaryHistos["MeanVDriftGlbSummary"] =
      ibooker.book2D("MeanVDriftGlbSummary", "mean VDrift average per sector", 12, 1., 13., 5, -2., 3.);
  allwheelHistos["allMeanVDrift"] =
      ibooker.book1D("VDriftMeanAllWheels", "mean VDrift for all chambers", 60, 0.0048, 0.006);

  ibooker.setCurrentFolder("DT/02-Segments/02-SigmaVDrift");

  summaryHistos["SigmaVDriftGlbSummary"] =
      ibooker.book2D("SigmaVDriftGlbSummary", "# of Chambers with good sigma VDrift", 12, 1., 13., 5, -2., 3.);
  allwheelHistos["allSigmaVDrift"] =
      ibooker.book1D("VDriftSigmaAllWheels", "sigma VDrift for all chambers", 30, 0., 0.0006);

  ibooker.setCurrentFolder("DT/02-Segments/03-MeanT0");

  summaryHistos["MeanT0GlbSummary"] =
      ibooker.book2D("MeanT0GlbSummary", "mean T0 average per sector", 12, 1., 13., 5, -2., 3.);
  allwheelHistos["allMeanT0"] = ibooker.book1D("T0MeanAllWheels", "mean T0 for all chambers", 100, -25., 25.);

  ibooker.setCurrentFolder("DT/02-Segments/03-SigmaT0");

  summaryHistos["SigmaT0GlbSummary"] =
      ibooker.book2D("SigmaT0GlbSummary", "# of Chambers with good sigma T0", 12, 1., 13., 5, -2., 3.);
  allwheelHistos["allSigmaT0"] = ibooker.book1D("T0SigmaAllWheels", "sigma T0 for all chambers", 50, 0, 25);

  for (int wh = -2; wh <= 2; wh++) {
    bookWheelHistos(ibooker, "MeanVDrift", "02-MeanVDrift", wh, 60, 0.0048, 0.006, true);
    bookWheelHistos(ibooker, "SigmaVDrift", "02-SigmaVDrift", wh, 30, 0., 0.0006);
    bookWheelHistos(ibooker, "MeanT0", "03-MeanT0", wh, 100, -25., 25., false, true);
    bookWheelHistos(ibooker, "SigmaT0", "03-SigmaT0", wh, 50, 0, 25, false, true);
  }

  for (int wheel = -2; wheel <= 2; wheel++) {
    for (int sec = 1; sec <= 14; sec++) {
      for (int stat = 1; stat <= 4; stat++) {
        if ((sec == 13 || sec == 14) && stat != 4)
          continue;

        // Get the ME produced by DTRunConditionVar Source
        MonitorElement* VDriftME = getChamberHistos(igetter, DTChamberId(wheel, stat, sec), "VDrift_FromSegm");
        MonitorElement* T0ME = getChamberHistos(igetter, DTChamberId(wheel, stat, sec), "T0_FromSegm");

        if (!VDriftME || !T0ME) {
          edm::LogWarning("DTRunConditionVarClient") << "ME not available" << std::endl;
          return;
        }

        // Get the means per chamber
        float vDriftMean = VDriftME->getMean();
        T0ME->setAxisRange(-15, 15);
        float t0Mean = T0ME->getMean();

        // Get the sigma per chamber
        float vDriftSigma = VDriftME->getRMS();
        float t0Sigma = T0ME->getRMS();

        if (VDriftME->getEntries() != 0) {
          allwheelHistos["allMeanVDrift"]->Fill(vDriftMean);
          allwheelHistos["allSigmaVDrift"]->Fill(vDriftSigma);

          (wheelHistos[wheel])["MeanVDrift"]->Fill(vDriftMean);
          (wheelHistos[wheel])["SigmaVDrift"]->Fill(vDriftSigma);
        }

        if (T0ME->getEntries() != 0) {
          allwheelHistos["allMeanT0"]->Fill(t0Mean);
          allwheelHistos["allSigmaT0"]->Fill(t0Sigma);

          (wheelRingHistos[wheel][stat])["MeanT0"]->Fill(t0Mean);
          (wheelRingHistos[wheel][stat])["SigmaT0"]->Fill(t0Sigma);
        }

        DTChamberId indexCh(wheel, stat, sec);

        float vDriftDev(0.), errvDriftDev(0.);
        percDevVDrift(indexCh, vDriftMean, vDriftSigma, vDriftDev, errvDriftDev);

        int sec_ = sec;
        if (sec == 13 || sec == 14)
          sec_ = (sec == 13) ? 4 : 10;

        float fillvDriftDev = max(min(vDriftDev, maxRangeVDrift), minRangeVDrift);
        float fillT0Mean = max(min(t0Mean, maxRangeT0), minRangeT0);

        float vDriftDevQ = varQuality(fabs(vDriftDev), maxGoodVDriftDev, minBadVDriftDev);
        float t0MeanQ = varQuality(fabs(t0Mean), maxGoodT0, minBadT0);

        float vDriftSigmQ = varQuality(vDriftSigma, maxGoodVDriftSigma, minBadVDriftSigma);
        float t0SigmQ = varQuality(t0Sigma, maxGoodT0Sigma, minBadT0Sigma);

        if (sec == 13 || sec == 14) {
          float binVDriftDev = (wheelHistos[wheel])["MeanVDriftSummary"]->getBinContent(sec_, stat);
          binVDriftDev = (fabs(binVDriftDev) > fabs(fillvDriftDev)) ? binVDriftDev : fillvDriftDev;
          (wheelHistos[wheel])["MeanVDriftSummary"]->setBinContent(sec_, stat, binVDriftDev);

          float binT0MeanVal = (wheelHistos[wheel])["MeanT0Summary"]->getBinContent(sec_, stat);
          binT0MeanVal = (fabs(binT0MeanVal) > fabs(fillT0Mean)) ? binT0MeanVal : fillT0Mean;
          (wheelHistos[wheel])["MeanT0Summary"]->setBinContent(sec_, stat, binT0MeanVal);

          float binVDriftSigmVal = (wheelHistos[wheel])["SigmaVDriftSummary"]->getBinContent(sec_, stat);
          binVDriftSigmVal = (binVDriftSigmVal > 0. && binVDriftSigmVal < vDriftSigmQ) ? binVDriftSigmVal : vDriftSigmQ;
          (wheelHistos[wheel])["SigmaVDriftSummary"]->setBinContent(sec_, stat, binVDriftSigmVal);

          float binT0SigmVal = (wheelHistos[wheel])["SigmaT0Summary"]->getBinContent(sec_, stat);
          binT0SigmVal = (binT0SigmVal > 0. && binT0SigmVal < t0SigmQ) ? binT0SigmVal : t0SigmQ;
          (wheelHistos[wheel])["SigmaT0Summary"]->setBinContent(sec_, stat, binT0SigmVal);

        } else {
          (wheelHistos[wheel])["MeanVDriftSummary"]->setBinContent(sec_, stat, fillvDriftDev);
          (wheelHistos[wheel])["MeanT0Summary"]->setBinContent(sec_, stat, fillT0Mean);
          (wheelHistos[wheel])["SigmaVDriftSummary"]->setBinContent(sec_, stat, vDriftSigmQ);
          (wheelHistos[wheel])["SigmaT0Summary"]->setBinContent(sec_, stat, t0SigmQ);
        }

        double weight = 1 / 4.;
        if ((sec_ == 4 || sec_ == 10) && stat == 4)
          weight = 1 / 8.;

        if (vDriftDevQ > 0.85 && vDriftSigmQ > 0.85) {
          glbVDriftSummary->Fill(sec_, wheel, weight);
          summaryHistos["MeanVDriftGlbSummary"]->Fill(sec_, wheel, weight);
          summaryHistos["SigmaVDriftGlbSummary"]->Fill(sec_, wheel, weight);

        } else {
          if (vDriftDevQ > 0.85 && vDriftSigmQ < 0.85) {
            summaryHistos["MeanVDriftGlbSummary"]->Fill(sec_, wheel, weight);
          }
          if (vDriftDevQ < 0.85 && vDriftSigmQ > 0.85) {
            summaryHistos["SigmaVDriftGlbSummary"]->Fill(sec_, wheel, weight);
          }
        }

        if (t0MeanQ > 0.85 && t0SigmQ > 0.85) {
          glbT0Summary->Fill(sec_, wheel, weight);
          summaryHistos["MeanT0GlbSummary"]->Fill(sec_, wheel, weight);
          summaryHistos["SigmaT0GlbSummary"]->Fill(sec_, wheel, weight);
        } else {
          if (t0MeanQ > 0.85 && t0SigmQ < 0.85) {
            summaryHistos["MeanT0GlbSummary"]->Fill(sec_, wheel, weight);
          }
          if (t0MeanQ < 0.85 && t0SigmQ > 0.85) {
            summaryHistos["SigmaT0GlbSummary"]->Fill(sec_, wheel, weight);
          }
        }

      }  // end loop on stations
    }  // end loop on sectors
  }  //end loop on wheels

  return;
}

float DTRunConditionVarClient::varQuality(float var, float maxGood, float minBad) {
  float qual(0);
  if (var <= maxGood) {
    qual = 1.;
  } else if (var > maxGood && var < minBad) {
    qual = 0.9;
  } else if (var >= minBad) {
    qual = 0.1;
  }

  return qual;
}

void DTRunConditionVarClient::percDevVDrift(
    DTChamberId indexCh, float meanVD, float sigmaVD, float& devVD, float& errdevVD) {
  DTSuperLayerId indexSLPhi1(indexCh, 1);
  DTSuperLayerId indexSLPhi2(indexCh, 3);

  float vDriftPhi1(0.), vDriftPhi2(0.);
  float ResPhi1(0.), ResPhi2(0.);
  if (readLegacyVDriftDB) {  // Legacy format
    int status1 = mTimeMap_->get(indexSLPhi1, vDriftPhi1, ResPhi1, DTVelocityUnits::cm_per_ns);
    int status2 = mTimeMap_->get(indexSLPhi2, vDriftPhi2, ResPhi2, DTVelocityUnits::cm_per_ns);

    if (status1 != 0 || status2 != 0) {
      DTSuperLayerId sl = (status1 != 0) ? indexSLPhi1 : indexSLPhi2;
      throw cms::Exception("DTRunConditionVarClient") << "Could not find vDrift entry in DB for" << sl << endl;
    }
  } else {
    vDriftPhi1 = vDriftMap_->get(DTWireId(indexSLPhi1.rawId()));
    vDriftPhi2 = vDriftMap_->get(DTWireId(indexSLPhi2.rawId()));
  }

  float vDriftMed = (vDriftPhi1 + vDriftPhi2) / 2.;

  devVD = (meanVD - vDriftMed) / vDriftMed;
  devVD = devVD < 1. ? devVD : 1.;

  errdevVD = sigmaVD / vDriftMed;

  return;
}

void DTRunConditionVarClient::bookWheelHistos(DQMStore::IBooker& ibooker,
                                              string histoType,
                                              string subfolder,
                                              int wh,
                                              int nbins,
                                              float min,
                                              float max,
                                              bool isVDCorr,
                                              bool makeRings) {
  stringstream wheel;
  wheel << wh;

  string folder = "DT/02-Segments/" + subfolder;

  ibooker.setCurrentFolder(folder);

  string histoName;
  string histoLabel;

  if (makeRings) {
    ibooker.setCurrentFolder(folder + "/Wheel" + wheel.str());
    for (int st = 1; st <= 4; st++) {
      stringstream station;
      station << st;

      histoName = histoType + "_W" + wheel.str() + "_MB" + station.str();
      histoLabel = histoType;

      (wheelRingHistos[wh][st])[histoType] = ibooker.book1D(histoName, histoLabel, nbins, min, max);
    }
  } else {
    histoName = histoType + "_W" + wheel.str();
    histoLabel = histoType;

    (wheelHistos[wh])[histoType] = ibooker.book1D(histoName, histoLabel, nbins, min, max);
  }

  ibooker.setCurrentFolder(folder);

  if (isVDCorr) {
    histoLabel = "Summary of corrections to VDrift DB values";
    histoName = "CorrTo" + histoType + "Summary_W" + wheel.str();
  } else {
    histoLabel = histoType + "Summary";
    histoName = histoType + "Summary_W" + wheel.str();
  }

  MonitorElement* me = ibooker.book2D(histoName, histoLabel, 12, 1, 13, 4, 1, 5);

  me->setBinLabel(1, "MB1", 2);
  me->setBinLabel(2, "MB2", 2);
  me->setBinLabel(3, "MB3", 2);
  me->setBinLabel(4, "MB4", 2);
  me->setAxisTitle("Sector", 1);

  (wheelHistos[wh])[histoType + "Summary"] = me;

  return;
}

DTRunConditionVarClient::MonitorElement* DTRunConditionVarClient::getChamberHistos(DQMStore::IGetter& igetter,
                                                                                   const DTChamberId& dtCh,
                                                                                   string histoType) {
  int wh = dtCh.wheel();
  int sc = dtCh.sector();
  int st = dtCh.station();
  stringstream wheel;
  wheel << wh;
  stringstream station;
  station << st;
  stringstream sector;
  sector << sc;

  string folder = "DT/02-Segments/Wheel" + wheel.str() + "/Sector" + sector.str() + "/Station" + station.str();
  string histoTag = "_W" + wheel.str() + "_Sec" + sector.str() + "_St" + station.str();
  string MEpath = folder + "/" + histoType + histoTag;

  igetter.setCurrentFolder(folder);

  LogTrace("DTDQM|DTMonitorModule|DTRunConditionVar") << "[DTRunConditionVar]: getting ME from " << folder << endl;

  MonitorElement* ME = igetter.get(MEpath);

  return ME;
}
