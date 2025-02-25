#ifndef SimG4Core_TmpSimTrack_H
#define SimG4Core_TmpSimTrack_H

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "SimG4Core/Notification/interface/TrackWithHistory.h"
#include <cmath>

class TmpSimTrack {
public:
  TmpSimTrack(int iid, int ipart, const math::XYZVectorD& ip, double ie)
      : id_(iid), ipart_(ipart), ip_(ip), ie_(ie), ivert_(-1), igenpart_(-1), parentID_(-1) {}

  TmpSimTrack(int iid, int ipart, const math::XYZVectorD& ip, double ie, int iv, int ig, const math::XYZVectorD& ipmom)
      : id_(iid), ipart_(ipart), ip_(ip), ie_(ie), ivert_(iv), igenpart_(ig), parentMomentum_(ipmom) {}

  TmpSimTrack(int iid,
              int ipart,
              const math::XYZVectorD& ip,
              double ie,
              int iv,
              int ig,
              const math::XYZVectorD& ipmom,
              const math::XYZVectorD& tkpos,
              const math::XYZTLorentzVectorD& tkmom)
      : id_(iid),
        ipart_(ipart),
        ip_(ip),
        ie_(ie),
        ivert_(iv),
        igenpart_(ig),
        parentMomentum_(ipmom),
        tkSurfacePosition_(tkpos),
        tkSurfaceMomentum_(tkmom) {}

  ~TmpSimTrack() = default;

  int id() const { return id_; }
  int part() const { return ipart_; }
  const math::XYZVectorD& momentum() const { return ip_; }
  double energy() const { return ie_; }
  int ivert() const { return ivert_; }
  int igenpart() const { return isPrimary_ ? igenpart_ : -1; }
  // parent momentum at interaction
  const math::XYZVectorD& parentMomentum() const { return parentMomentum_; }
  // Information at level of tracker surface
  const math::XYZVectorD& trackerSurfacePosition() const { return tkSurfacePosition_; }
  const math::XYZTLorentzVectorD& trackerSurfaceMomentum() const { return tkSurfaceMomentum_; }
  // parent track ID (only stored if parent momentum at interaction
  // is stored, else = -1)
  int parentID() const { return parentID_; }

  void copyCrossedBoundaryVars(const TrackWithHistory* track) {
    if (track->crossedBoundary()) {
      crossedBoundary_ = track->crossedBoundary();
      idAtBoundary_ = track->getIDAtBoundary();
      positionAtBoundary_ = track->getPositionAtBoundary();
      momentumAtBoundary_ = track->getMomentumAtBoundary();
    }
  }
  bool crossedBoundary() const { return crossedBoundary_; }
  const math::XYZTLorentzVectorF& getPositionAtBoundary() const { return positionAtBoundary_; }
  const math::XYZTLorentzVectorF& getMomentumAtBoundary() const { return momentumAtBoundary_; }
  int getIDAtBoundary() const { return idAtBoundary_; }
  bool isFromBackScattering() const { return isFromBackScattering_; }
  void setFromBackScattering() { isFromBackScattering_ = true; }
  void setGenParticleID(int i) { igenpart_ = i; }
  bool isPrimary() const { return isPrimary_; }
  void setIsPrimary() { isPrimary_ = true; }
  int getPrimaryID() const { return igenpart_; }

private:
  int id_;
  int ipart_;
  math::XYZVectorD ip_;
  double ie_;
  int ivert_;
  int igenpart_;
  int parentID_;
  math::XYZVectorD parentMomentum_{math::XYZVectorD(0., 0., 0.)};
  math::XYZVectorD tkSurfacePosition_{math::XYZVectorD(0., 0., 0.)};
  math::XYZTLorentzVectorD tkSurfaceMomentum_{math::XYZTLorentzVectorD(0., 0., 0., 0.)};
  bool isFromBackScattering_{false};
  bool crossedBoundary_{false};
  bool isPrimary_{false};
  int idAtBoundary_{-1};
  math::XYZTLorentzVectorF positionAtBoundary_{math::XYZTLorentzVectorF(0.f, 0.f, 0.f, 0.f)};
  math::XYZTLorentzVectorF momentumAtBoundary_{math::XYZTLorentzVectorF(0.f, 0.f, 0.f, 0.f)};
};

#endif
