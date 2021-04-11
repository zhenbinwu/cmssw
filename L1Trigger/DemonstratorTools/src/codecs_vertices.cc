
#include "L1Trigger/DemonstratorTools/interface/codecs/vertices.h"

namespace l1t::demo::codecs {

  ap_uint<64> encodeVertex(const l1t::VertexWord& v) { return v.vertexWord(); }

  // Encodes vertex collection onto 1 output link
  std::array<l1t::demo::BoardData::Channel, 1> encodeVertices(const edm::View<l1t::VertexWord>& vertices) {
    std::vector<ap_uint<64>> vertexWords;

    for (const auto& vertex : vertices)
      vertexWords.push_back(encodeVertex(vertex));

    std::array<l1t::demo::BoardData::Channel, 1> linkData;

    for (size_t i = 0; i < linkData.size(); i++) {
      // Pad vertex vectors -> full packet length (10 frames = 10 vertices)
      vertexWords.resize(10, 0);
      linkData.at(i).resize(vertexWords.size(), {0});

      for (size_t j = 0; (j < vertexWords.size()); j++) {
        linkData.at(i).at(j).data = vertexWords.at(j)(63, 0);
      }
    }

    return linkData;
  }

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
