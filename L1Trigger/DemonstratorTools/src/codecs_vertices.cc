
#include "L1Trigger/DemonstratorTools/interface/codecs/vertices.h"

namespace l1t::demo::codecs {

  std::vector<l1t::Vertex> decodeVertices(const std::vector<l1t::demo::Frame>& frames) {
    std::vector<l1t::Vertex> vertices;

    for (const auto& x : frames) {
      if (not(x.data & 1))
        break;

      // TODO: Replace with complete implementation
      Vertex vertex(ap_fixed<15, 6>(ap_uint<64>(x.data)(15, 1)), {});
      vertices.push_back(vertex);
    }

    return vertices;
  }

}  // namespace l1t::demo::codecs