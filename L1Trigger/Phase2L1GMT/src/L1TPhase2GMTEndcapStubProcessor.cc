#include "L1Trigger/Phase2L1GMT/interface/L1TPhase2GMTEndcapStubProcessor.h"
#include "L1Trigger/L1TMuon/interface/MuonTriggerPrimitive.h"
#include <cmath>
#include <iostream> 
#include <string> 
#include <sstream> 

L1TPhase2GMTEndcapStubProcessor::L1TPhase2GMTEndcapStubProcessor():
  minBX_(-3),
  maxBX_(3)
{
 
} 



L1TPhase2GMTEndcapStubProcessor::L1TPhase2GMTEndcapStubProcessor(const edm::ParameterSet& iConfig):
  minBX_(iConfig.getParameter<int>("minBX")),
  maxBX_(iConfig.getParameter<int>("maxBX")),
  coord1LSB_(iConfig.getParameter<double>("coord1LSB")),
  coord2LSB_(iConfig.getParameter<double>("coord2LSB")),
  eta1LSB_(iConfig.getParameter<double>("eta1LSB")),
  eta2LSB_(iConfig.getParameter<double>("eta2LSB")),
  etaMatch_(iConfig.getParameter<double>("etaMatch")),
  phiMatch_(iConfig.getParameter<double>("phiMatch")),
  verbose_(iConfig.getParameter<unsigned int>("verbose"))
{

} 



L1TPhase2GMTEndcapStubProcessor::~L1TPhase2GMTEndcapStubProcessor() {}




L1TPhase2GMTStub
L1TPhase2GMTEndcapStubProcessor::buildCSCOnlyStub(const CSCDetId& detid,const CSCCorrelatedLCTDigi& digi, const L1TMuon::GeometryTranslator *translator) {

  int endcap = detid.endcap();
  int station = detid.station();
  int chamber =detid.chamber();
  int ring = detid.ring();

  L1TMuon::TriggerPrimitive primitive(detid,digi);

  const GlobalPoint& gp = translator->getGlobalPoint(primitive);


  int phi = int(gp.phi().value()/coord1LSB_);
  int eta1 = int(gp.eta()/eta1LSB_);


  int wheel=0;
  int sign = endcap==1 ? -1 : 1;

  if (ring==3)
    wheel = sign*3;
  if (ring==2)
    wheel = sign*4;
  if (ring==1)
    wheel = sign*5;

  int sector = fabs(chamber);

  int bx=digi.getBX()-8;
  int quality=1;

  uint tfLayer=0;


  L1TPhase2GMTStub stub(wheel,sector,station,tfLayer,phi,0,0,
			bx,quality,eta1,0,1,0); 

  stub.setOfflineQuantities(gp.phi().value(),0.0,gp.eta(),0.0);
  return stub;

  }


L1TPhase2GMTStub
L1TPhase2GMTEndcapStubProcessor::buildRPCOnlyStub(const RPCDetId& detid,const RPCDigi& digi, const L1TMuon::GeometryTranslator *translator) {


  L1TMuon::TriggerPrimitive primitive(detid,digi);
  const GlobalPoint& gp = translator->getGlobalPoint(primitive);


  int phi2 = int(gp.phi().value()/coord2LSB_);
  int eta2 = int(gp.eta()/eta2LSB_);

  int wheel=-(6-detid.ring())*detid.region();
  int sector =(detid.sector()-1)*6+detid.subsector();
  int station=detid.station();
  bool tag = detid.trIndex();
  int bx=digi.bx();
  int quality=2;
  uint tfLayer=0;

  L1TPhase2GMTStub stub(wheel,sector,station,tfLayer,0,phi2,tag,
			bx,quality,0,eta2,2,0); 
  stub.setOfflineQuantities(0.0,gp.phi().value(),0.0,gp.eta());
  return stub;

  }


L1TPhase2GMTStubCollection 
L1TPhase2GMTEndcapStubProcessor::combineStubs(const L1TPhase2GMTStubCollection& cscStubs,const L1TPhase2GMTStubCollection& rpcStubs) {
  L1TPhase2GMTStubCollection out;
  L1TPhase2GMTStubCollection usedRPC;
  L1TPhase2GMTStubCollection cleanedRPC;

  for (const auto & csc : cscStubs) {
    int nRPC=0;
    float phiF=0.0;
    float etaF=0.0;
    int phi=0;
    int eta=0;
    for (const auto & rpc : rpcStubs) {
      if (csc.etaRegion()!=rpc.etaRegion() || csc.depthRegion()!=rpc.depthRegion())
	continue;
      if (fabs(deltaPhi(csc.offline_coord1(),rpc.offline_coord2()))<phiMatch_ &&
	  fabs(csc.offline_eta1()-rpc.offline_eta2())<etaMatch_ && csc.bxNum()==rpc.bxNum()) {
	phiF+=rpc.offline_coord2();
	etaF+=rpc.offline_eta2();
	phi+=rpc.coord2();
	eta+=rpc.eta2();
	nRPC++;
	usedRPC.push_back(rpc);
      }

    }
    
    int finalRPCPhi=0;
    int finalRPCEta=0;
    double offline_finalRPCPhi=0;
    double offline_finalRPCEta=0;
    if (nRPC!=0) {
      finalRPCPhi=phi/nRPC;
      finalRPCEta=eta/nRPC;
      offline_finalRPCPhi=phiF/nRPC;
      offline_finalRPCEta=etaF/nRPC;
      
    }


    //make the fancy stub
    L1TPhase2GMTStub stub(csc.etaRegion(),csc.phiRegion(),csc.depthRegion(),0,csc.coord1(),finalRPCPhi,0,
			  csc.bxNum(),nRPC|(1<<3),csc.eta1(),finalRPCEta,2,0); 
    stub.setOfflineQuantities(csc.offline_coord1(),offline_finalRPCPhi,csc.offline_eta1(),offline_finalRPCEta);
    out.push_back(stub);
  }

    //clean the RPC from the used ones
    
  for (const auto & rpc : rpcStubs) {
    bool keep=true;
    for (const auto & rpc2 : usedRPC) {
      if (rpc==rpc2) {
	keep=false;
	break;
      }
    }
    if (keep)
      cleanedRPC.push_back(rpc);
  }

  while(cleanedRPC.size()>0) {
    L1TPhase2GMTStubCollection freeRPC;

    int nRPC=1;
    float phiF=cleanedRPC[0].offline_coord2();
    float etaF=cleanedRPC[0].offline_eta2();
    int phi=cleanedRPC[0].coord2();
    int eta=cleanedRPC[0].eta2();

    for (unsigned i=1;i<cleanedRPC.size();++i) {
      if (fabs(deltaPhi(cleanedRPC[0].offline_coord2(),cleanedRPC[i].offline_coord2()))<phiMatch_ &&
	  cleanedRPC[0].etaRegion()==cleanedRPC[i].etaRegion() && cleanedRPC[0].depthRegion()==cleanedRPC[i].depthRegion() &&
	  fabs(cleanedRPC[0].offline_eta2()-cleanedRPC[i].offline_eta2())<etaMatch_ &&cleanedRPC[0].bxNum()==cleanedRPC[i].bxNum()) {
	  phiF+=cleanedRPC[i].offline_coord2();
	  etaF+=cleanedRPC[i].offline_eta2();
	  phi+=cleanedRPC[i].coord2();
	  eta+=cleanedRPC[i].eta2();
	  nRPC++;
      } 
      else {
	freeRPC.push_back(cleanedRPC[i]);
      }
    }
    L1TPhase2GMTStub stub(cleanedRPC[0].etaRegion(),cleanedRPC[0].phiRegion(),cleanedRPC[0].depthRegion(),0,0,phi/nRPC,0,
			  cleanedRPC[0].bxNum(),nRPC,0,eta/nRPC,2,0); 
    stub.setOfflineQuantities(0.0,phiF/nRPC,0.0,etaF/nRPC);
    out.push_back(stub);
    cleanedRPC=freeRPC;
  };
  return out;
}




L1TPhase2GMTStubCollection 
L1TPhase2GMTEndcapStubProcessor::makeStubs(const MuonDigiCollection<CSCDetId,CSCCorrelatedLCTDigi>& csc,const MuonDigiCollection<RPCDetId,RPCDigi>& cleaned,const L1TMuon::GeometryTranslator *t,const edm::EventSetup& iSetup) {
  L1TPhase2GMTStubCollection cscStubs;
   auto chamber = csc.begin();
   auto chend  = csc.end();
   for( ; chamber != chend; ++chamber ) {
     auto digi = (*chamber).second.first;
     auto dend = (*chamber).second.second;    
     for( ; digi != dend; ++digi ) {
       L1TPhase2GMTStub stub = buildCSCOnlyStub((*chamber).first,*digi,t);
       if (stub.bxNum()>=minBX_&&stub.bxNum()<=maxBX_)
	 cscStubs.push_back(stub);
     }
   }



  L1TPhase2GMTStubCollection rpcStubs;

  //  RPCHitCleaner cleaner(rpc);
  //cleaner.run(iSetup);    
  //  RPCDigiCollection cleaned = rpc;

   auto rpcchamber = cleaned.begin();
   auto rpcchend  = cleaned.end();
   for( ; rpcchamber != rpcchend; ++rpcchamber ) {
       if((*rpcchamber).first.region()==0)
	 continue;
     auto digi = (*rpcchamber).second.first;
     auto dend = (*rpcchamber).second.second;    
     for( ; digi != dend; ++digi ) {
       L1TPhase2GMTStub stub = buildRPCOnlyStub((*rpcchamber).first,*digi,t);
       if (stub.bxNum()>=minBX_&&stub.bxNum()<=maxBX_)
	 rpcStubs.push_back(stub);
     }
   }


  L1TPhase2GMTStubCollection combinedStubs =combineStubs(cscStubs,rpcStubs);


   if (verbose_) {
     printf("CSC Stubs\n");
     for (const auto& stub : cscStubs) 
       printf("CSC Stub bx=%d etaRegion=%d phiRegion=%d depthRegion=%d  coord1=%f,%d coord2=%f,%d eta1=%f,%d eta2=%f,%d quality=%d etaQuality=%d\n",
	      stub.bxNum(),stub.etaRegion(),stub.phiRegion(),stub.depthRegion(),stub.offline_coord1(),stub.coord1(),stub.offline_coord2(),stub.coord2(),stub.offline_eta1(),stub.eta1(),stub.offline_eta2(),stub.eta2(),stub.quality(),stub.etaQuality());
     printf("RPC Stubs\n");
     for (const auto& stub : rpcStubs) 
       printf("RPC Stub bx=%d etaRegion=%d phiRegion=%d depthRegion=%d  coord1=%f,%d coord2=%f,%d eta1=%f,%d eta2=%f,%d quality=%d etaQuality=%d\n",
	      stub.bxNum(),stub.etaRegion(),stub.phiRegion(),stub.depthRegion(),stub.offline_coord1(),stub.coord1(),stub.offline_coord2(),stub.coord2(),stub.offline_eta1(),stub.eta1(),stub.offline_eta2(),stub.eta2(),stub.quality(),stub.etaQuality());
     for (const auto& stub : combinedStubs) 
       printf("Combined Stub bx=%d etaRegion=%d phiRegion=%d depthRegion=%d  coord1=%f,%d coord2=%f,%d eta1=%f,%d eta2=%f,%d quality=%d etaQuality=%d\n",
	      stub.bxNum(),stub.etaRegion(),stub.phiRegion(),stub.depthRegion(),stub.offline_coord1(),stub.coord1(),stub.offline_coord2(),stub.coord2(),stub.offline_eta1(),stub.eta1(),stub.offline_eta2(),stub.eta2(),stub.quality(),stub.etaQuality());

   }
  
   return combinedStubs;
}

