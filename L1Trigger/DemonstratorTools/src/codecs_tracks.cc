
#include "L1Trigger/DemonstratorTools/interface/codecs/tracks.h"

namespace l1t::demo::codecs {

  ap_uint<96> encodeTrack(const TTTrack_TrackWord& t) {
    // FIXME: Update to actually encode track parameters (just hardcoding values
    //        right now to verify file-writing and packing of 96-bit words).

    ap_uint<96> word(0);
    word.set(0);   // valid
    word.set(1);   // pt
    word.set(16);  // phi
    word.set(28);  // tanLambda
    word.set(44);  // z0
    word.set(56);  // d0
    word.set(69);  // chi2 (r-phi)
    word.set(73);  // chi2 (r-z)
    word.set(77);  // bend chi2
    word.set(80);  // hit pattern
    word.set(87);  // MVA track quality
    word.set(90);  // Other MVAs

    return word;
  }

  // Encodes track collection onto 18 output links (2x9 eta-phi sectors; first 9 negative eta)
  std::array<std::vector<l1t::demo::Frame>, 18> encodeTracks(const edm::View<TTTrack<Ref_Phase2TrackerDigi_>>& tracks) {
    std::array<std::vector<ap_uint<96>>, 18> trackWords;

    for (const auto& track : tracks)
      trackWords.at((track.eta() >= 0 ? 9 : 0) + track.phiSector()).push_back(encodeTrack(track));

    std::array<std::vector<l1t::demo::Frame>, 18> linkData;

    for (size_t i = 0; i < linkData.size(); i++) {
      // Pad track vectors -> full packet length (156 frames = 104 tracks)
      trackWords.at(i).resize(104, 0);
      linkData.at(i).resize(156, {0});

      for (size_t j = 0; (j < trackWords.size()); j += 2) {
        linkData.at(i).at(3 * j / 2).data = ap_uint<64>(trackWords.at(i).at(j)(63, 0));
        linkData.at(i).at(3 * j / 2 + 1).data =
            ap_uint<64>((ap_uint<32>(trackWords.at(i).at(j + 1)(31, 0)), ap_uint<32>(trackWords.at(i).at(j)(95, 64))));
        linkData.at(i).at(3 * j / 2 + 2).data = ap_uint<64>(trackWords.at(i).at(j + 1)(95, 32));
      }
    }

    return linkData;
  }

}  // namespace l1t::demo::codecs