#ifndef L1Trigger_L1TTrackMatch_Cordic_HH
#define L1Trigger_L1TTrackMatch_Cordic_HH

#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuAlgo.h"

using namespace L1TkEtMissEmuAlgo;

class Cordic {
public:
  Cordic();
  Cordic(iPhi aPhiScale, int aMagnitudeBits, const int aSteps, bool debug, bool writeLUTs);

  EtMiss to_polar(iEt x, iEt y) const;

private:
  const int mPhiScale;
  const int mMagnitudeScale;
  const int mMagnitudeBits;
  const int cordic_steps;

  const bool debug;

  std::vector<iPhi> atan_LUT;
  std::vector<iEt> mag_renormalization_LUT;
};

#endif