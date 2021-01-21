#include "L1Trigger/Phase2L1GMT/interface/L1TPhase2GMTBarrelStubProcessor.h"
#include <cmath>
#include <iostream> 
#include <string> 
#include <sstream> 

L1TPhase2GMTBarrelStubProcessor::L1TPhase2GMTBarrelStubProcessor():
  minPhiQuality_(0),
  minBX_(-3),
  maxBX_(3)
{
 
} 



L1TPhase2GMTBarrelStubProcessor::L1TPhase2GMTBarrelStubProcessor(const edm::ParameterSet& iConfig):
  minPhiQuality_(iConfig.getParameter<int>("minPhiQuality")),
  minBX_(iConfig.getParameter<int>("minBX")),
  maxBX_(iConfig.getParameter<int>("maxBX")),
  eta1_(iConfig.getParameter<std::vector<int> >("eta_1")),
  eta2_(iConfig.getParameter<std::vector<int> >("eta_2")),
  eta3_(iConfig.getParameter<std::vector<int> >("eta_3")),
  coarseEta1_(iConfig.getParameter<std::vector<int> >("coarseEta_1")),
  coarseEta2_(iConfig.getParameter<std::vector<int> >("coarseEta_2")),
  coarseEta3_(iConfig.getParameter<std::vector<int> >("coarseEta_3")),
  coarseEta4_(iConfig.getParameter<std::vector<int> >("coarseEta_4")),
  verbose_(iConfig.getParameter<int>("verbose")),
  phiLSB_(iConfig.getParameter<double>("phiLSB"))

{

} 



L1TPhase2GMTBarrelStubProcessor::~L1TPhase2GMTBarrelStubProcessor() {}








L1TPhase2GMTStub 
L1TPhase2GMTBarrelStubProcessor::buildStub(const L1Phase2MuDTPhDigi& phiS,const L1MuDTChambThDigi* etaS) {
  
  L1TPhase2GMTStub stub = buildStubNoEta(phiS);


  //Now full eta
  int qeta1=0;
  int qeta2=0;
  int eta1=255;
  int eta2=255; 


  bool hasEta=false;
  for (uint i=0;i<7;++i) {
    if (etaS->position(i)==0)
      continue;
    if (!hasEta) {
      hasEta=true;
      eta1=calculateEta(i,etaS->whNum(),etaS->scNum(),etaS->stNum());
      if (etaS->quality(i)==1)
	qeta1=2;
      else
	qeta1=1;
    }
    else {
      eta2=calculateEta(i,etaS->whNum(),etaS->scNum(),etaS->stNum());
      if (etaS->quality(i)==1)
	qeta2=2;
      else
	qeta2=1;
    }
  }



  if (qeta2>0) {//both stubs->average
    stub.setEta(eta1,eta2,qeta1+qeta2);
  }
  else if (qeta1>0) {//Good single stub
    stub.setEta(eta1,0,qeta1);
  }
 


  return stub;

}





L1TPhase2GMTStub
L1TPhase2GMTBarrelStubProcessor::buildStubNoEta(const L1Phase2MuDTPhDigi& phiS) {
  int wheel = phiS.whNum();
  int abswheel = fabs(phiS.whNum());
  int sign  = wheel>0 ? 1: -1;
  int sector = phiS.scNum();
  int station = phiS.stNum();
  double globalPhi = (sector*30)+phiS.phi()*30./2048.;
  if (globalPhi<-180)
    globalPhi+=360;
  if (globalPhi>180)
    globalPhi-=360;
  globalPhi = globalPhi*M_PI/180.;
  int phi  = int(globalPhi/phiLSB_);
  int phiB = phiS.phiBend();
  uint tag = phiS.index();
  int bx=phiS.bxNum();
  int quality=phiS.quality();
  uint tfLayer=0;
  int eta=-255;
  if (station==1) {
    eta=-coarseEta1_[abswheel];
  }
  else if (station==2) {
    eta=-coarseEta2_[abswheel];
  }
  else if (station==3) {
    eta=-coarseEta3_[abswheel];
  }
  else if (station==4) {
    eta=-coarseEta4_[abswheel];
  }
  else {
    eta=-255;
  }
  //Now full eta

  eta=eta*sign;
  L1TPhase2GMTStub stub(wheel,sector,station,tfLayer,phi,phiB,tag,
			bx,quality,eta,0,0,0);
  stub.setOfflineQuantities(globalPhi,float(phiB),0.0,0.0);
  return stub;
}





L1TPhase2GMTStubCollection 
L1TPhase2GMTBarrelStubProcessor::makeStubs(const L1Phase2MuDTPhContainer* phiContainer,const L1MuDTChambThContainer* etaContainer) {

  L1TPhase2GMTStubCollection out;
  for (int bx=minBX_;bx<=maxBX_;bx++) {
    for (int wheel=-2;wheel<=2;wheel++) {
      for (int sector=0;sector<12;sector++) {
	for (int station=1;station<5;station++) {

	  bool hasEta=false;
	  L1MuDTChambThDigi const*  tseta   = etaContainer->chThetaSegm(wheel,station,sector,bx);
	  if (tseta) {
	    hasEta=true;
	  }
	  for (const auto& phiDigi : *phiContainer->getContainer()) {
	    if (phiDigi.bxNum()!= bx ||phiDigi.whNum()!=wheel || phiDigi.scNum()!=sector || phiDigi.stNum()!=station)
	      continue;
	    if (hasEta) {
	      out.push_back(buildStub(phiDigi,tseta));
	    }
	    else {
	      out.push_back(buildStubNoEta(phiDigi));
	    }
	  }
	}
      }
    }
  }

  return out;
}



int L1TPhase2GMTBarrelStubProcessor::calculateEta(uint i, int wheel,uint sector,uint station) {
  int eta=0;
  if (wheel>0) {
	eta=7*wheel+3-i;
      }
  else if (wheel<0) {
	eta=7*wheel+i-3;
  }
  else {
    if (sector==0 || sector==3 ||sector==4 ||sector==7 ||sector==8 ||sector==11)
      eta=i-3;
    else
      eta=3-i;
  }


  if (station==1)
    eta=-eta1_[eta+17];
  else if (station==2)
    eta=-eta2_[eta+17];
  else 
    eta=-eta3_[eta+17];

  return eta;



}



