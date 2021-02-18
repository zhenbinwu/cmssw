
#ifndef L1Trigger_DemonstratorTools_BoardDataReader_h
#define L1Trigger_DemonstratorTools_BoardDataReader_h

#include <map>
#include <string>
#include <vector>

#include "L1Trigger/DemonstratorTools/interface/BoardData.h"
#include "L1Trigger/DemonstratorTools/interface/ChannelSpec.h"
#include "L1Trigger/DemonstratorTools/interface/FileFormat.h"

namespace l1t::demo {

  // Reads I/O buffer files created from hardware/firmware tests, verifying that
  // received packets conform to the declared structure, separating out each
  // event (accounting for different TM periods of specific links and of the
  // data-processor itself), and transparently switching to data from new buffer
  // files as needed
  class BoardDataReader {
  public:
    BoardDataReader(FileFormat,
                    const std::vector<std::string>&,
                    const size_t framesPerBX,
                    const size_t tmux,
                    const size_t emptyFramesAtStart,
                    const std::map<size_t, ChannelSpec>&);

    BoardData getNextEvent();

  private:
    FileFormat fileFormat_;

    std::vector<std::string> fileNames_;

    size_t framesPerBX_;

    size_t boardTMUX_;

    size_t emptyFramesAtStart_;

    std::map<size_t, ChannelSpec> channelSpecs_;

    std::vector<BoardData> events_;

    std::vector<BoardData>::const_iterator eventIt_;
  };

}  // namespace l1t::demo

#endif