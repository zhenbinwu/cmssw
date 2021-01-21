#ifndef L1TPHASE2GMTBARRELSTUBPROCESSOR
#define L1TPHASE2GMTBARRELSTUBPROCESSOR

#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/L1TMuonPhase2/interface/L1TPhase2GMTStub.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/L1TObjects/interface/L1TMuonBarrelParams.h"
#include "CondFormats/DataRecord/interface/L1TMuonBarrelParamsRcd.h"




class L1TPhase2GMTBarrelStubProcessor {
 public:
  L1TPhase2GMTBarrelStubProcessor();
  L1TPhase2GMTBarrelStubProcessor(const edm::ParameterSet&);
  
  ~L1TPhase2GMTBarrelStubProcessor();


  L1TPhase2GMTStubCollection makeStubs(const L1Phase2MuDTPhContainer*,const L1MuDTChambThContainer*);
  
 private:
  L1TPhase2GMTStub buildStub(const L1Phase2MuDTPhDigi&,const L1MuDTChambThDigi*);
  L1TPhase2GMTStub buildStubNoEta(const L1Phase2MuDTPhDigi&);

  int calculateEta(uint, int,uint,uint);  
  int minPhiQuality_;
  int minBX_;
  int maxBX_;
  std::vector<int> eta1_;
  std::vector<int> eta2_;
  std::vector<int> eta3_;
  std::vector<int> coarseEta1_;
  std::vector<int> coarseEta2_;
  std::vector<int> coarseEta3_;
  std::vector<int> coarseEta4_;
  int verbose_;
  double phiLSB_;

};


#endif
