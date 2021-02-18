
#ifndef L1Trigger_DemonstratorTools_utilities_h
#define L1Trigger_DemonstratorTools_utilities_h

#include "L1Trigger/DemonstratorTools/interface/BoardData.h"
#include "L1Trigger/DemonstratorTools/interface/FileFormat.h"

namespace l1t::demo {

  BoardData read(const std::string& filePath, const FileFormat format);

  void writeToFile(const BoardData& data, const std::string& filePath, const FileFormat format);

}  // namespace l1t::demo

#endif