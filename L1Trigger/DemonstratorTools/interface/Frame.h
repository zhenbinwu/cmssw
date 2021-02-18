
#ifndef L1Trigger_DemonstratorTools_Frame_h
#define L1Trigger_DemonstratorTools_Frame_h

#include <cstdint>

#include "ap_int.h"

namespace l1t::demo {

  struct Frame {
    Frame() = default;

    Frame(const uint64_t);

    Frame(const ap_uint<64>&);

    uint64_t data{0};
    bool valid{false};
    bool strobe{true};
  };

}  // namespace l1t::demo

#endif
