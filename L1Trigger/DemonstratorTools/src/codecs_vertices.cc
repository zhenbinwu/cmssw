
#include "L1Trigger/DemonstratorTools/interface/codecs/vertices.h"

namespace l1t::demo::codecs {

  std::vector<l1t::VertexWord> decodeVertices(const l1t::demo::BoardData::Channel& frames) {
    std::vector<l1t::VertexWord> vertices;

    for (const auto& x : frames) {
      if (not x.data.test(VertexWord::kValidLSB))
        break;

      VertexWord v(VertexWord::vtxvalid_t(1),
                   VertexWord::vtxz0_t(x.data(VertexWord::kZ0MSB, VertexWord::kZ0LSB)),
                   VertexWord::vtxmultiplicity_t(x.data(VertexWord::kNTrackInPVMSB, VertexWord::kNTrackInPVLSB)),
                   VertexWord::vtxsumpt_t(x.data(VertexWord::kSumPtMSB, VertexWord::kSumPtLSB)),
                   VertexWord::vtxquality_t(x.data(VertexWord::kQualityMSB, VertexWord::kQualityLSB)),
                   VertexWord::vtxinversemult_t(x.data(VertexWord::kNTrackOutPVMSB, VertexWord::kNTrackOutPVLSB)),
                   VertexWord::vtxunassigned_t(x.data(VertexWord::kUnassignedMSB, VertexWord::kUnassignedLSB)));
      vertices.push_back(v);
    }

    return vertices;
  }

}  // namespace l1t::demo::codecs