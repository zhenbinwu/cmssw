#include "DataFormats/L1TMuonPhase2/interface/L1TPhase2GMTStub.h"

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
L1TPhase2GMTStub::L1TPhase2GMTStub() :
  etaRegion_(0),phiRegion_(0),depthRegion_(0),coord1_(0), coord2_(0),id_(0), quality_(-1), bxNum_(17),
  eta1_(0),eta2_(0),etaQuality_(-1),type_(0) {}

L1TPhase2GMTStub::L1TPhase2GMTStub(int etaRegion,int phiRegion,int depthRegion,uint tfLayer,int coord1,int coord2,int id,int bx,int quality,int eta1,int eta2,int etaQuality,int type):
  etaRegion_(etaRegion),
  phiRegion_(phiRegion),
  depthRegion_(depthRegion),
  tfLayer_(tfLayer),
  coord1_(coord1),
  coord2_(coord2),
  id_(id),
  quality_(quality),
  bxNum_(bx),
  eta1_(eta1),
  eta2_(eta2),
  etaQuality_(etaQuality),
  type_(type)
{

}

L1TPhase2GMTStub::~L1TPhase2GMTStub() {}


bool L1TPhase2GMTStub::operator==(const L1TPhase2GMTStub& id) const {

  if ( etaRegion_             != id.etaRegion_ )               return false;
  if ( phiRegion_             != id.phiRegion_ )               return false;
  if ( depthRegion_           != id.depthRegion_)               return false;
  if ( id_                    != id.id_ )                      return false;
  if ( coord1_                != id.coord1_ )                  return false;
  if ( coord2_                != id.coord2_ )                  return false;
  if ( quality_               != id.quality_ )                 return false;
  if ( bxNum_                 != id.bxNum_ )                   return false;
  if ( eta1_                  != id.eta1_ )                    return false;
  if ( eta2_                  != id.eta2_ )                    return false;
  if ( etaQuality_            != id.etaQuality_ )              return false;
  if ( type_                  != id.type_ )                    return false;
  return true;
}

//
// output stream operator for phi track segments
//
ostream& operator<<(ostream& s, const L1TPhase2GMTStub& id) {

  s.setf(ios::right,ios::adjustfield);
  s << "BX: "              << setw(5) << id.bxNum_  << " "
  << "etaregion:"          << setw(5) << id.etaRegion_  << " "
  << "phiRegion: "         << setw(5) << id.phiRegion_  << " "
  << "depthRegion: "       << setw(5) << id.depthRegion_  << " "
  << "stub ID: "           << setw(5) << id.id_  << " "
  << "coord1: "            << setw(5) << id.coord1_  << " "
  << "coord2: "            << setw(4) << id.coord2_ << " "
  << "quality: "           << setw(4) << id.quality_ << " "
  << "eta1:"               << setw(4) << id.eta1_ << " "
  << "eta2:"               << setw(4) << id.eta2_ << " "
  << "qeta:"               << setw(4) << id.etaQuality_ << " "
  << "type:"               << setw(4) <<id.type_;
  return s;

}
