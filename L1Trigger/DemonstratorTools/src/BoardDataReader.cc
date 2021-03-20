#include "L1Trigger/DemonstratorTools/interface/BoardDataReader.h"

#include <fstream>

#include "L1Trigger/DemonstratorTools/interface/Frame.h"
#include "L1Trigger/DemonstratorTools/interface/utilities.h"

namespace l1t::demo {

  BoardDataReader::BoardDataReader(FileFormat format,
                                   const std::vector<std::string>& fileNames,
                                   const size_t framesPerBX,
                                   const size_t tmux,
                                   const size_t emptyFramesAtStart,
                                   const std::map<size_t, ChannelSpec>& channelSpecs)
      : fileFormat_(format),
        fileNames_(fileNames),
        framesPerBX_(framesPerBX),
        boardTMUX_(tmux),
        emptyFramesAtStart_(emptyFramesAtStart),
        channelSpecs_(channelSpecs),
        events_() {
    // TODO: Move much of this to separate function, and only read files on demand
    for (const auto& path : fileNames_) {
      BoardData boardData(read(path, fileFormat_));

      // 1) Verify that all expected channels are present
      for (const auto& [i, chanSpec] : channelSpecs_) {
        if (not boardData.has(i))
          throw std::runtime_error("Channel " + std::to_string(i) + " was declared but is missing from file '" + path +
                                   "'");
      }

      // 2) Verify that packet structure is as expected
      for (const auto& [i, chanSpec] : channelSpecs_) {
        const auto& chanData = boardData.at(i);

        const size_t framesBeforeFirstPacket(emptyFramesAtStart_ + chanSpec.tmuxIndex * framesPerBX_ + chanSpec.offset);
        const size_t eventLength(chanSpec.tmux * framesPerBX_);
        const size_t packetLength(eventLength - chanSpec.interpacketGap);

        for (size_t j = 0; j < framesBeforeFirstPacket; j++) {
          if (chanData.at(j).valid)
            throw std::runtime_error("Frame " + std::to_string(j) + " on channel " + std::to_string(i) +
                                     " is valid, but first " + std::to_string(framesBeforeFirstPacket) +
                                     "frames should be invalid");
        }

        for (size_t j = framesBeforeFirstPacket; j < chanData.size(); j++) {
          bool expectValid(((j - framesBeforeFirstPacket) % eventLength) < packetLength);

          if (expectValid) {
            if (not chanData.at(j).valid)
              throw std::runtime_error("Frame " + std::to_string(j) + " on channel " + std::to_string(i) +
                                       " is invalid, but expected valid frame");
          } else if (chanData.at(j).valid)
            throw std::runtime_error("Frame " + std::to_string(j) + " on channel " + std::to_string(i) +
                                     " is valid, but expected invalid frame");
        }
      }

      // 3) Extract the data for each event
      bool eventIncomplete(false);
      for (size_t eventIndex = 0;; eventIndex++) {
        BoardData eventData;

        for (const auto& [i, chanSpec] : channelSpecs_) {
          const auto& chanData = boardData.at(i);

          // Skip this channel if its data is not part of this event
          if ((chanSpec.tmuxIndex / boardTMUX_) != ((eventIndex / chanSpec.tmux) / boardTMUX_))
            continue;

          // Extract the frames for this event
          const size_t framesBeforeEvent((chanSpec.tmuxIndex + eventIndex * boardTMUX_) * framesPerBX_ +
                                         emptyFramesAtStart_ + chanSpec.offset);
          const size_t packetLength(chanSpec.tmux * framesPerBX_ - chanSpec.interpacketGap);

          if (chanData.size() < (framesBeforeEvent + chanSpec.tmux * framesPerBX_)) {
            eventIncomplete = true;
            break;
          }

          BoardData::Channel chanEventData(packetLength);
          for (size_t j = 0; j < packetLength; j++)
            chanEventData.at(j) = chanData.at(framesBeforeEvent + j);
          eventData.add(i, chanEventData);
        }

        if (eventIncomplete)
          break;

        events_.push_back(eventData);
      }
    }

    eventIt_ = events_.begin();
  }

  BoardData BoardDataReader::getNextEvent() {
    if (eventIt_ == events_.end())
      throw std::runtime_error("Board data reader ran out of events");

    return *(eventIt_++);
  }

}  // namespace l1t::demo