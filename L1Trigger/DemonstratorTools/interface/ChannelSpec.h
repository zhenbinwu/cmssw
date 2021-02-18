#ifndef L1Trigger_DemonstratorTools_ChannelSpec_h
#define L1Trigger_DemonstratorTools_ChannelSpec_h

#include <stddef.h>

namespace l1t::demo {

  struct ChannelSpec {
  public:
    // Time multiplexing period of data on link
    size_t tmux;
    // Index of events handled on link (divided by board's TMUX)
    size_t tmuxIndex;
    // Number of invalid frames between packets (i.e. following each event)
    size_t interpacketGap;
    // Number of invalid frames before first valid packet, beyond tmuxIndex * framesPerBX
    size_t offset{0};
  };

}  // namespace l1t::demo

#endif