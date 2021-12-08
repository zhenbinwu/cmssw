#include "l2egencoder_ref.h"

using namespace l1ct;

#ifdef CMSSW_GIT_HASH

#include "FWCore/ParameterSet/interface/ParameterSet.h"

l1ct::L2EgEncoderEmulator::L2EgEncoderEmulator(const edm::ParameterSet &pset)
    : L2EgEncoderEmulator(pset.getParameter<uint32_t>("n64BitWords")) {}

#endif

void L2EgEncoderEmulator::toFirmware(const std::vector<ap_uint<64>>& encoded_in, ap_uint<64> encoded_fw[]) const {
  for(unsigned int i = 0; i < nEncodedWords_; i++) {
    encoded_fw[i] = (i < encoded_in.size()) ? encoded_in[i] : ap_uint<64>(0);
  }
} 

std::vector<ap_uint<64>> L2EgEncoderEmulator::encodeLayer2EgObjs(unsigned int nObj, 
  const std::vector<EGIsoObjEmu>& photons, 
  const std::vector<EGIsoEleObjEmu>& electrons) const {
    std::vector<ap_uint<64>> ret;
    
    auto encoded_photons = encodeLayer2(photons);
    encoded_photons.resize(nObj, {0});
    auto encoded_eles = encodeLayer2(electrons);
    encoded_eles.resize(nObj, {0});

    encodeLayer2To64bits(encoded_photons, ret);
    encodeLayer2To64bits(encoded_eles, ret);

    return ret;
}

void L2EgEncoderEmulator::encodeLayer2To64bits(const std::vector<ap_uint<96>>& packed96, std::vector<ap_uint<64>>& packed64) const {
  for(unsigned int i = 0; i < packed96.size(); i+=2) {

    packed64.push_back(packed96[i](63, 0));
    packed64.push_back((ap_uint<32>(packed96[i+1](95, 64)), ap_uint<32>(packed96[i](95, 64))));
    packed64.push_back(packed96[i+1](63, 0));

    // std::cout << "obj [" << i << "]: " << std::hex << packed96[i] << std::endl;
    // std::cout << "obj [" << i+1 << "]: " << std::hex << packed96[i+1] << std::endl;
    // std::cout << "frame [" << std::dec << packed64.size()-3 << "]" << std::hex << packed64[packed64.size()-3] << std::endl;
    // std::cout << "frame [" << std::dec << packed64.size()-2 << "]" << std::hex << packed64[packed64.size()-2] << std::endl;
    // std::cout << "frame [" << std::dec << packed64.size()-1 << "]" << std::hex << packed64[packed64.size()-1] << std::endl;
    
  }
}
