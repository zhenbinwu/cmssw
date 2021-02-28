#include <memory>
#include <cmath>

#include "L1Trigger/L1TTrackMatch/interface/Cordic.h"


using namespace L1TkEtMissEmuAlgo;

Cordic::Cordic(iPhi aPhiScale, int aMagnitudeBits, int aSteps) :
mPhiScale(aPhiScale),mMagnitudeScale(1 << aMagnitudeBits), mMagnitudeBits(aMagnitudeBits), cordic_steps(aSteps){

    atan_LUT.reserve(cordic_steps);
    mag_renormalization_LUT.reserve(cordic_steps);

    for (int i = 0; i< cordic_steps; i++){
        atan_LUT.push_back(iPhi(round(mPhiScale*atan(pow(2,-i))/(2*M_PI))));
    }

    float val = 1.0;
    for (int j =0; j<cordic_steps;j++){
        val = val / (pow(1 + pow(4,-j),0.5));
        mag_renormalization_LUT.push_back(iEt(round(mMagnitudeScale * val)));
      }  
  
}

EtMiss Cordic::to_polar(iEt x,iEt y) const {

    iEt new_x = 0;
    iEt new_y = 0;

    iglobPhi phi = 0;
    iglobPhi new_phi = 0;
    bool sign = false;

    EtMiss ret_etmiss;

    if (x >=0 &&  y >= 0){
            phi = 0;
            sign = true;
            x = x;
            y = y;
        }
        else if( x < 0 &&  y >= 0){
            phi = phiScale >> 1;
            sign = false;
            x = -x;
            y = - y;
        }
        else if( x < 0 &&  y < 0){
            phi = phiScale >> 1;
            sign = true;
            x = -x;
            y = - y;
        }
        else{
            phi = phiScale;
            sign = false;
            x = x;
            y = - y;
        }
        
        for (int step = 0; step<cordic_steps; step++){
            if (y < 0){
              new_x = x - (y >> step);
              new_y = y + (x >> step);
            } else {
              new_x = x + (y >> step);
              new_y = y - (x >> step); 
            }

            if ((y < 0) == sign){
                new_phi = phi - atan_LUT[step];
            } else{
                new_phi = phi + atan_LUT[step];
            }      

            x = new_x;
            y = new_y;
            phi = new_phi;

        }
        
        ret_etmiss.Et = x*mag_renormalization_LUT[cordic_steps-1] >> mMagnitudeBits;
        ret_etmiss.Phi = iglobPhi (phi);
        return ret_etmiss;
      
}


