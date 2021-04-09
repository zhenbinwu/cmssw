
#ifndef L1Trigger_DemonstratorTools_codecs_vertices_h
#define L1Trigger_DemonstratorTools_codecs_vertices_h

#include <array>
#include <vector>

#include "ap_int.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/L1Trigger/interface/Vertex.h"
#include "DataFormats/L1Trigger/interface/VertexWord.h"
#include "L1Trigger/DemonstratorTools/interface/BoardData.h"

namespace l1t::demo::codecs {

  ap_uint<64> encodeVertex(const l1t::VertexWord& v);

  // Encodes vertex collection onto 1 'logical' output link
  std::array<l1t::demo::BoardData::Channel, 1> encodeVertices(const edm::View<l1t::VertexWord>&);

  std::vector<l1t::VertexWord> decodeVertices(const l1t::demo::BoardData::Channel&);

}  // namespace l1t::demo::codecs

#endif
