#ifndef L2EGENCODER_REF_H
#define L2EGENCODER_REF_H

#ifdef CMSSW_GIT_HASH
#include "../dataformats/layer1_emulator.h"
#include "../dataformats/egamma.h"
#include "L1Trigger/Phase2L1ParticleFlow/src/dbgPrintf.h"

#else
#include "../../dataformats/layer1_emulator.h"
#include "../../dataformats/egamma.h"
#include "../../utils/dbgPrintf.h"

#endif


namespace edm {
  class ParameterSet;
}


namespace l1ct {

  struct L2EgEncoderEmulator {
    public:

      L2EgEncoderEmulator(unsigned int nEncodedWords) : nEncodedWords_(nEncodedWords) {};
      
      L2EgEncoderEmulator(const edm::ParameterSet &iConfig);

      void toFirmware(const std::vector<ap_uint<64>>& encoded_in, ap_uint<64> encoded_fw[]) const;


      
      std::vector<ap_uint<64>> encodeLayer2EgObjs(unsigned int nObj, 
        const std::vector<EGIsoObjEmu>& photons, 
        const std::vector<EGIsoEleObjEmu>& electrons) const;

    private:

      template<class T>
      ap_uint<96> encodeLayer2(const T& egiso) const {
        ap_uint<96> ret = 0;
        // FIXME; should be packed in GT format
        ret(T::BITWIDTH, 0) = egiso.pack();
        return ret;
      }


      template<class T>
      std::vector<ap_uint<96>> encodeLayer2(const std::vector<T>& egisos) const {
        std::vector<ap_uint<96>> ret;
        for(const auto&egiso: egisos) {
          ret.push_back(encodeLayer2(egiso));
        }
        return ret;
      }

      void encodeLayer2To64bits(const std::vector<ap_uint<96>>& packed96, std::vector<ap_uint<64>>& packed64) const;
      
      unsigned int nEncodedWords_;
  };

}
#endif
