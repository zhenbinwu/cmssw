
#ifndef L1Trigger_DemonstratorTools_codecs_vertices_h
#define L1Trigger_DemonstratorTools_codecs_vertices_h

#include <vector>

#include "DataFormats/L1Trigger/interface/VertexWord.h"

#include "L1Trigger/DemonstratorTools/interface/BoardData.h"

namespace l1t::demo::codecs {

  std::vector<l1t::VertexWord> decodeVertices(const l1t::demo::BoardData::Channel&);

}

#endif