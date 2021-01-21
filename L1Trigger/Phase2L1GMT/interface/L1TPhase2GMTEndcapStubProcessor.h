
#ifndef L1TPHASE2GMTENDCAPSTUBPROCESSOR
#define L1TPHASE2GMTENDCAPSTUBPROCESSOR

#include "DataFormats/L1TMuonPhase2/interface/L1TPhase2GMTStub.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "L1Trigger/L1TMuon/interface/GeometryTranslator.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "L1Trigger/L1TMuon/interface/MuonTriggerPrimitive.h"
#include "L1Trigger/L1TTwinMux/interface/RPCHitCleaner.h"

class L1TPhase2GMTEndcapStubProcessor {


 public:
  L1TPhase2GMTEndcapStubProcessor();
  L1TPhase2GMTEndcapStubProcessor(const edm::ParameterSet&);
    ~L1TPhase2GMTEndcapStubProcessor();



    L1TPhase2GMTStubCollection makeStubs(const MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>& csc,const MuonDigiCollection<RPCDetId,RPCDigi>& rpc,const L1TMuon::GeometryTranslator *t,const edm::EventSetup& iSetup);

 private:
    L1TPhase2GMTStub buildCSCOnlyStub(const CSCDetId& , const CSCCorrelatedLCTDigi&,const L1TMuon::GeometryTranslator*);
    L1TPhase2GMTStub buildRPCOnlyStub(const RPCDetId& , const RPCDigi&,const L1TMuon::GeometryTranslator*);
    L1TPhase2GMTStubCollection combineStubs(const L1TPhase2GMTStubCollection&,const L1TPhase2GMTStubCollection&);
    

  int minBX_;
  int maxBX_;
  double coord1LSB_;
  double coord2LSB_;
  double eta1LSB_;
  double eta2LSB_;
  double etaMatch_;
  double phiMatch_;
  bool verbose_;


};


#endif
