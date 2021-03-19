#include "L1Trigger/DemonstratorTools/interface/BoardDataWriter.h"

#include <fstream>

#include "L1Trigger/DemonstratorTools/interface/Frame.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"

namespace l1t::demo {

  BoardDataWriter::BoardDataWriter(FileFormat format,
                                   const std::string& path,
                                   const size_t framesPerBX,
                                   const size_t tmux,
                                   const size_t maxFramesPerFile,
                                   const std::map<size_t, ChannelSpec>& channelSpecs)
      : fileFormat_(format),
        filePathGen_([=](const size_t i) { return path + "_" + std::to_string(i) + ".txt"; }),
        framesPerBX_(framesPerBX),
        boardTMUX_(tmux),
        maxFramesPerFile_(maxFramesPerFile),
        maxEventsPerFile_(maxFramesPerFile_),
        eventIndex_(0),
        channelSpecs_(channelSpecs) {
    for (const auto& [i, x] : channelSpecs_) {
      boardData_.add(i);
      maxEventsPerFile_ = std::min(
          maxEventsPerFile_, (x.tmux / boardTMUX_) * size_t((maxFramesPerFile_ - x.tmuxIndex * framesPerBX_ - x.offset) / (x.tmux * framesPerBX_)));
    }

    resetBoardData();
  }

  void BoardDataWriter::addEvent(const BoardData& eventData) {
    for (const auto& [i, channelData] : eventData) {
      if (channelSpecs_.count(i) == 0)
        throw std::runtime_error("Event data for link " + std::to_string(i) +
                                 " was given to BoardDataFileWriter, but its structure was not defined");

      const auto& chanSpec = channelSpecs_.at(i);
      if (channelData.size() > (channelSpecs_.at(i).tmux * framesPerBX_ - chanSpec.interpacketGap))
        throw std::runtime_error("Event data for link " + std::to_string(i) + " (TMUX " +
                                 std::to_string(chanSpec.tmux) + ", " + std::to_string(chanSpec.interpacketGap) +
                                 " frames between packets) is too long (" + std::to_string(channelData.size()) +
                                 " frames)");

      boardData_.at(i).insert(boardData_.at(i).end(), channelData.begin(), channelData.end());
      boardData_.at(i).insert(
          boardData_.at(i).end(), channelSpecs_.at(i).tmux * framesPerBX_ - channelData.size(), Frame());
    }

    eventIndex_++;

    if ((eventIndex_ % maxEventsPerFile_) == 0) {
      // Pad any channels that aren't full with invalid frames
      for (auto& x : boardData_)
        x.second.resize(maxFramesPerFile_);

      // Write board data object to file
      const std::string filePath = filePathGen_(fileNames_.size());
      writeToFile(boardData_, filePath, fileFormat_);
      fileNames_.push_back(filePath);

      // Clear board data to be ready for next event
      resetBoardData();
    }
  }

  void BoardDataWriter::resetBoardData() {
    for (auto& x : boardData_)
      x.second.clear();

    for (const auto& [i, spec] : channelSpecs_)
      boardData_.at(i).resize(spec.tmuxIndex * framesPerBX_ + spec.offset);
  }

}  // namespace l1t::demo