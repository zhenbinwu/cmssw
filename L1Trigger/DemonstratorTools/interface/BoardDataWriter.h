
#ifndef L1Trigger_DemonstratorTools_BoardDataWriter_h
#define L1Trigger_DemonstratorTools_BoardDataWriter_h

#include <functional>

#include "L1Trigger/DemonstratorTools/interface/BoardData.h"
#include "L1Trigger/DemonstratorTools/interface/ChannelSpec.h"
#include "L1Trigger/DemonstratorTools/interface/FileFormat.h"

namespace l1t::demo {

  // Writes I/O buffer files created from hardware/firmware tests, ensuring that
  // the data conforms to the declared packet structure (by inserting invalid
  // frames automatically), concatenating data from each event (accounting for
  // different TM periods of specific links, and of the data-processor itself),
  // and transparently switching to new output files when the data would overrun
  // the length of the board's I/O buffers
  class BoardDataWriter {
  public:
    BoardDataWriter(FileFormat,
                    const std::string& filePath,
                    const size_t framesPerBX,
                    const size_t tmux,
                    const size_t maxFramesPerFile,
                    const std::map<size_t, ChannelSpec>&);

    void addEvent(const BoardData& data);

  private:
    void resetBoardData();

    FileFormat fileFormat_;

    std::function<std::string(const size_t)> filePathGen_;

    std::vector<std::string> fileNames_;

    size_t framesPerBX_;

    size_t boardTMUX_;

    size_t maxFramesPerFile_;

    size_t maxEventsPerFile_;

    size_t eventIndex_;

    BoardData boardData_;

    // map of link channel index -> TMUX period, TMUX index, interpacket-gap & offset
    std::map<size_t, ChannelSpec> channelSpecs_;
  };

}  // namespace l1t::demo

#endif