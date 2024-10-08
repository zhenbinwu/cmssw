#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "CondFormats/HcalObjects/interface/HcalRecoParams.h"
#include "CondFormats/DataRecord/interface/HcalRecoParamsRcd.h"

#include "CondFormats/HcalObjects/interface/HcalPFCuts.h"
#include "CondFormats/DataRecord/interface/HcalPFCutsRcd.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloTopology/interface/HcalTopology.h"
#include "Geometry/HcalTowerAlgo/interface/HcalGeometry.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalChannelProperties.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalChannelPropertiesAuxRecord.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalChannelPropertiesRecord.h"

class HcalChannelPropertiesEP : public edm::ESProducer {
public:
  typedef std::unique_ptr<HcalRecoParams> ReturnType1;
  typedef std::unique_ptr<HcalChannelPropertiesVec> ReturnType2;
  typedef std::unique_ptr<HcalPFCuts> ReturnType3;

  inline HcalChannelPropertiesEP(const edm::ParameterSet&) {
    auto cc1 = setWhatProduced(this, &HcalChannelPropertiesEP::produce1);
    topoToken_ = cc1.consumes();
    paramsToken_ = cc1.consumes();

    auto cc2 = setWhatProduced(this, &HcalChannelPropertiesEP::produce2);
    edm::ESInputTag qTag("", "withTopo");
    condToken_ = cc2.consumes();
    myParamsToken_ = cc2.consumes();
    sevToken_ = cc2.consumes();
    qualToken_ = cc2.consumes(qTag);
    geomToken_ = cc2.consumes();

    edm::es::Label productLabel("withTopo");
    auto cc3 = setWhatProduced(this, &HcalChannelPropertiesEP::produce3, productLabel);
    topoToken3_ = cc3.consumes();
    pfcutsToken_ = cc3.consumes();
  }

  inline ~HcalChannelPropertiesEP() override {}

  ReturnType1 produce1(const HcalChannelPropertiesAuxRecord& rcd) {
    const HcalTopology& htopo = rcd.get(topoToken_);
    const HcalRecoParams& params = rcd.get(paramsToken_);

    ReturnType1 prod = std::make_unique<HcalRecoParams>(params);
    prod->setTopo(&htopo);
    return prod;
  }

  ReturnType2 produce2(const HcalChannelPropertiesRecord& rcd) {
    // There appears to be no easy way to trace the internal
    // dependencies of HcalDbService. So, rebuild the product
    // every time anything changes in the parent records.
    // This means that we are sometimes going to rebuild the
    // whole table on the lumi block boundaries instead of
    // just updating the list of bad channels.
    //
    // Retrieve various event setup records and data products
    const HcalDbRecord& dbRecord = rcd.getRecord<HcalDbRecord>();
    const HcalDbService& cond = dbRecord.get(condToken_);
    const HcalRecoParams& params = rcd.get(myParamsToken_);
    const HcalSeverityLevelComputer& severity = rcd.get(sevToken_);
    const HcalChannelQuality& qual = dbRecord.get(qualToken_);
    const CaloGeometry& geom = rcd.get(geomToken_);

    // HcalTopology is taken from "params" created by the "produce1" method
    const HcalTopology& htopo(*params.topo());

    // Build the product
    ReturnType2 prod = std::make_unique<HcalChannelPropertiesVec>(htopo.ncells());
    std::array<HcalPipelinePedestalAndGain, 4> pedsAndGains;
    const HcalSubdetector subdetectors[3] = {HcalBarrel, HcalEndcap, HcalForward};

    for (HcalSubdetector subd : subdetectors) {
      const HcalGeometry* hcalGeom = static_cast<const HcalGeometry*>(geom.getSubdetectorGeometry(DetId::Hcal, subd));
      const std::vector<DetId>& ids = hcalGeom->getValidDetIds(DetId::Hcal, subd);

      for (const auto cell : ids) {
        const auto rawId = cell.rawId();

        // ADC decoding tools, etc
        const HcalRecoParam* param_ts = params.getValues(rawId);
        const HcalQIECoder* channelCoder = cond.getHcalCoder(cell);
        const HcalQIEShape* shape = cond.getHcalShape(channelCoder);
        const HcalSiPMParameter* siPMParameter = cond.getHcalSiPMParameter(cell);

        // Pedestals and gains
        const HcalCalibrations& calib = cond.getHcalCalibrations(cell);
        const HcalCalibrationWidths& calibWidth = cond.getHcalCalibrationWidths(cell);
        for (int capid = 0; capid < 4; ++capid) {
          pedsAndGains[capid] = HcalPipelinePedestalAndGain(calib.pedestal(capid),
                                                            calibWidth.pedestal(capid),
                                                            calib.effpedestal(capid),
                                                            calibWidth.effpedestal(capid),
                                                            calib.respcorrgain(capid),
                                                            calibWidth.gain(capid));
        }

        // Channel quality
        const HcalChannelStatus* digistatus = qual.getValues(rawId);
        const bool taggedBadByDb = severity.dropChannel(digistatus->getValue());

        // Fill the table entry
        const unsigned linearId = htopo.detId2denseId(cell);
        prod->at(linearId) =
            HcalChannelProperties(&calib, param_ts, channelCoder, shape, siPMParameter, pedsAndGains, taggedBadByDb);
      }
    }

    return prod;
  }

  ReturnType3 produce3(const HcalPFCutsRcd& rcd) {
    const HcalTopology& htopo = rcd.get(topoToken3_);
    const HcalPFCuts& cuts = rcd.get(pfcutsToken_);

    ReturnType3 prod = std::make_unique<HcalPFCuts>(cuts);
    prod->setTopo(&htopo);
    return prod;
  }

  HcalChannelPropertiesEP() = delete;
  HcalChannelPropertiesEP(const HcalChannelPropertiesEP&) = delete;
  HcalChannelPropertiesEP& operator=(const HcalChannelPropertiesEP&) = delete;

private:
  edm::ESGetToken<HcalDbService, HcalDbRecord> condToken_;
  edm::ESGetToken<HcalTopology, HcalRecNumberingRecord> topoToken_, topoToken3_;
  edm::ESGetToken<HcalRecoParams, HcalRecoParamsRcd> paramsToken_;
  edm::ESGetToken<HcalSeverityLevelComputer, HcalSeverityLevelComputerRcd> sevToken_;
  edm::ESGetToken<HcalChannelQuality, HcalChannelQualityRcd> qualToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;
  edm::ESGetToken<HcalRecoParams, HcalChannelPropertiesAuxRecord> myParamsToken_;
  edm::ESGetToken<HcalPFCuts, HcalPFCutsRcd> pfcutsToken_;
};

DEFINE_FWK_EVENTSETUP_MODULE(HcalChannelPropertiesEP);
