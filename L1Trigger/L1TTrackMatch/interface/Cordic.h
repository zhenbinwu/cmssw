#ifndef L1Trigger_L1TTrackMatch_Cordic_HH
#define L1Trigger_L1TTrackMatch_Cordic_HH

#include "L1Trigger/L1TTrackMatch/interface/L1TkEtMissEmuAlgo.h"

using namespace L1TkEtMissEmuAlgo;

class Cordic {
public:
  Cordic();
  Cordic(int aPhiScale, int aMagnitudeBits, const int aSteps, bool debug, bool writeLUTs);

  EtMiss toPolar(iEt x, iEt y) const;

private:
  //Scale for Phi calculation to maintain precision
  const int mPhiScale;
  //Scale for Magnitude calculation
  const int mMagnitudeScale;
  //Bit width for internal magnitude
  const int mMagnitudeBits;
  //Number of cordic iterations
  const int cordicSteps;

  const bool debug;

  //To calculate atan
  std::vector<iMETphi> atanLUT;
  //To normalise final magnitude
  std::vector<iEt> magNormalisationLUT;
};

#endif