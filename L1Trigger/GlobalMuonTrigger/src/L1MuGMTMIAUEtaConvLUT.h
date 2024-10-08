//-------------------------------------------------
//
/** \class L1MuGMTMIAUEtaConvLUT
 *
 *   MIAUEtaConv look-up table
 *          
 *   this class was automatically generated by 
 *     L1MuGMTLUT::MakeSubClass()  
*/
//
//   Author :
//   H. Sakulin            HEPHY Vienna
//
//   Migrated to CMSSW:
//   I. Mikulec
//
//--------------------------------------------------
#ifndef L1TriggerGlobalMuonTrigger_L1MuGMTMIAUEtaConvLUT_h
#define L1TriggerGlobalMuonTrigger_L1MuGMTMIAUEtaConvLUT_h

//---------------
// C++ Headers --
//---------------

//----------------------
// Base Class Headers --
//----------------------
#include "L1Trigger/GlobalMuonTrigger/src/L1MuGMTLUT.h"

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

//              ---------------------
//              -- Class Interface --
//              ---------------------

class L1MuGMTMIAUEtaConvLUT : public L1MuGMTLUT {
public:
  enum { MIP_DT, MIP_BRPC, ISO_DT, ISO_BRPC, MIP_CSC, MIP_FRPC, ISO_CSC, ISO_FRPC };

  /// constuctor using function-lookup
  L1MuGMTMIAUEtaConvLUT()
      : L1MuGMTLUT("MIAUEtaConv",
                   "MIP_DT MIP_BRPC ISO_DT ISO_BRPC MIP_CSC MIP_FRPC ISO_CSC ISO_FRPC",
                   "eta_in(6)",
                   "eta_out(4)",
                   6,
                   true) {
    InitParameters();
  };

  /// destructor
  ~L1MuGMTMIAUEtaConvLUT() override {}

  /// specific lookup function for eta_out
  unsigned SpecificLookup_eta_out(int idx, unsigned eta_in) const {
    std::vector<unsigned> addr(1);
    addr[0] = eta_in;
    return Lookup(idx, addr)[0];
  };

  /// specific lookup function for entire output field
  unsigned SpecificLookup(int idx, unsigned eta_in) const {
    std::vector<unsigned> addr(1);
    addr[0] = eta_in;
    return LookupPacked(idx, addr);
  };

  /// access to lookup function with packed input and output

  unsigned LookupFunctionPacked(int idx, unsigned address) const override {
    std::vector<unsigned> addr = u2vec(address, m_Inputs);
    return TheLookupFunction(idx, addr[0]);
  };

private:
  /// Initialize scales, configuration parameters, alignment constants, ...
  void InitParameters();

  /// The lookup function - here the functionality of the LUT is implemented
  unsigned TheLookupFunction(int idx, unsigned eta_in) const;

  /// Private data members (LUT parameters);
};
#endif
