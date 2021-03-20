#include "L1Trigger/DemonstratorTools/interface/utilities.h"

#include <fstream>
#include <unordered_map>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

#include "L1Trigger/DemonstratorTools/interface/BoardData.h"

namespace {
  std::string searchForID(std::istream& file) {
    std::string line, id;

    while (getline(file, line)) {
      boost::trim(line);
      if (line.empty())
        continue;
      if (line[0] == '#')
        continue;

      if (line.rfind("Board ", 0) != std::string::npos)
        return line.substr(6);
      else
        throw std::logic_error("Found unexpected line found when searching for board ID: \"" + line + "\"");
    }
    throw std::logic_error("Board ID not found!");
  }

  std::vector<std::string> searchAndTokenize(std::istream& file, const std::string& linePrefix) {
    std::string line;

    while (getline(file, line)) {
      boost::trim(line);
      if (line.empty())
        continue;
      if (line[0] == '#')
        continue;

      if (line.rfind(linePrefix, 0) != std::string::npos) {
        std::vector<std::string> tokens;
        std::string lineContents(line.substr(linePrefix.size()));
        // Trim the line
        boost::trim(lineContents);
        // Split the line into tokens
        boost::split(tokens, lineContents, boost::is_any_of(" \t"), boost::token_compress_on);
        return tokens;
      } else
        throw std::logic_error("Found unexpected line found when searching for \"" + linePrefix + "\": \"" + line +
                               "\"");
    }
    throw std::logic_error("Couldn't find any line starting with \"" + linePrefix + "\"");
  }

  std::vector<size_t> searchForLinks(std::istream& file) {
    searchAndTokenize(file, "Quad/Chan :");
    const auto tokens = searchAndTokenize(file, "Link :");
    std::vector<size_t> links;
    std::transform(
        tokens.begin(), tokens.end(), std::back_inserter(links), boost::lexical_cast<size_t, const std::string&>);
    return links;
  }

  l1t::demo::Frame convertStringToFrame(const std::string& token) {
    static const boost::regex frameRegex("([01]s)?([01]v)([0-9a-fA-F]{16})");

    boost::smatch what;
    if (!boost::regex_match(token, what, frameRegex))
      throw std::logic_error("Token '" + token + "' doesn't match the valid format");

    l1t::demo::Frame value;
    // Import strobe if the strobe group is matched
    if (what[1].matched) {
      value.strobe = (what[1] == "1s");
    }

    value.valid = (what[2] == "1v");
    value.data = std::stoull(what[3].str(), nullptr, 16);

    return value;
  }

  std::vector<std::vector<l1t::demo::Frame>> readDataRows(std::istream& file) {
    std::string line;
    std::vector<std::vector<l1t::demo::Frame>> data;

    while (file.good() and getline(file, line)) {
      boost::trim(line);

      if (line.empty() or line[0] == '#')
        continue;

      std::ostringstream prefixStream;
      prefixStream << "Frame ";
      prefixStream << std::setw(4) << std::setfill('0') << data.size();
      prefixStream << " :";

      const std::string prefix(prefixStream.str());
      if (line.rfind(prefix, 0) != std::string::npos) {
        std::vector<std::string> tokens;
        std::string tmp(line.substr(prefix.size()));
        boost::trim(tmp);
        boost::split(tokens, tmp, boost::is_any_of(" \t"), boost::token_compress_on);

        std::vector<l1t::demo::Frame> row;
        std::transform(tokens.begin(), tokens.end(), std::back_inserter(row), convertStringToFrame);

        data.push_back(row);
      } else
        throw std::logic_error("Found unexpected line found when searching for \"" + prefix + "\": \"" + line + "\"");
    }

    return data;
  }
}  // namespace

// APx sideband encoding
//   Short-term, simulation only:
//     0 -> Valid
//     1 -> EOF
//   Planned (from ~ May 2021)
//     0 -> Valid
//     1 -> SOF (Start Of Frame)
//     2 -> FFO (First Frame of Orbit)
//     3 -> EOF (End Of Frame)
//     4 -> FERR (Frame Error)
//     5 -> RSV1
//     6 -> RSV2
//     7 -> RSV3

namespace l1t::demo {

  FileFormat parseFileFormat(const std::string& s) {
    static const std::unordered_map<std::string, FileFormat> kFormatStringMap({{"EMP", FileFormat::EMP},
                                                                               {"emp", FileFormat::EMP},
                                                                               {"APx", FileFormat::APx},
                                                                               {"apx", FileFormat::APx},
                                                                               {"X20", FileFormat::X20},
                                                                               {"x20", FileFormat::X20}});

    const auto it = kFormatStringMap.find(s);
    if (it == kFormatStringMap.end())
      throw std::runtime_error("Could not convert '" + s + "' to FileFormat enum value");

    return it->second;
  }

  BoardData readAPxFile(std::istream&, const FileFormat);

  BoardData readEMPFile(std::istream&, const FileFormat);

  BoardData readX20File(std::istream&, const FileFormat);

  BoardData read(const std::string& filePath, const FileFormat format) {
    std::ifstream file(filePath);

    if (not file.is_open())
      throw std::runtime_error("Could not open file '" + filePath + "'");

    return read(file, format);
  }

  BoardData read(std::istream& file, const FileFormat format) {
    if (format == FileFormat::APx)
      return readAPxFile(file, format);
    else if (format == FileFormat::EMP)
      return readEMPFile(file, format);
    else
      return readX20File(file, format);
  }

  BoardData readAPxFile(std::istream& file, const FileFormat format) {
    std::cout << "WARNING: Reading APx file format not yet implemented. Will be done ASAP." << std::endl;
    return BoardData();
  }

  BoardData readEMPFile(std::istream& file, const FileFormat format) {
    std::string id = searchForID(file);
    BoardData boardData(id);

    std::vector<size_t> channels = searchForLinks(file);
    std::vector<std::vector<Frame>> dataRows = readDataRows(file);

    std::vector<std::vector<Frame>> dataColumns(channels.size(), std::vector<Frame>(dataRows.size()));

    for (size_t i = 0; i < channels.size(); i++)
      for (size_t j = 0; j < dataRows.size(); j++)
        dataColumns.at(i).at(j) = dataRows.at(j).at(i);

    for (size_t i = 0; i < channels.size(); i++)
      boardData.add(channels.at(i), dataColumns.at(i));

    return boardData;
  }

  BoardData readX20File(std::istream& file, const FileFormat format) {
    throw std::runtime_error("Reading X20 file format not yet implemented. Will be done ASAP.");
  }

  void writeAPxFile(const BoardData&, std::ostream&, const FileFormat);

  void writeEMPFile(const BoardData&, std::ostream&, const FileFormat);

  void writeX20File(const BoardData&, std::ostream&, const FileFormat);

  void write(const BoardData& data, const std::string& filePath, const FileFormat format) {
    // Open file
    std::cout << "Writing board data (" << std::distance(data.begin(), data.end()) << " channels, "
              << data.begin()->second.size() << " frames) to file '" << filePath << "' (format: " << format << ")"
              << std::endl;
    std::ofstream file(filePath);

    if (not file.is_open())
      throw std::runtime_error("Could not open file '" + filePath + "'");

    write(data, file, format);
  }

  void write(const BoardData& data, std::ostream& file, const FileFormat format) {
    // Check that number of frames is same for every channel
    const auto firstChannel = data.begin();

    for (const auto& [i, channelData] : data) {
      if (channelData.size() != firstChannel->second.size())
        throw std::runtime_error("Cannot write board data to file - channels do not all have the same length (" +
                                 std::to_string(channelData.size()) + " words on channel " + std::to_string(i) +
                                 ", but " + std::to_string(firstChannel->second.size()) + " words on channel " +
                                 std::to_string(firstChannel->first) + ")");
    }

    // Call relevant write function
    switch (format) {
      case FileFormat::APx:
        writeAPxFile(data, file, format);
        return;
      case FileFormat::EMP:
        writeEMPFile(data, file, format);
        return;
      case FileFormat::X20:
        writeX20File(data, file, format);
        return;
    }
  }

  void writeAPxFile(const BoardData& data, std::ostream& file, const FileFormat format) {
    file << std::setfill('0');

    file << "#Sideband ON" << std::endl;

    // Channel header
    file << "#LinkLabel";
    for (const auto& [i, channelData] : data)
      file << "                LINK_" << std::setw(2) << i << "    ";
    file << std::endl;

    file << "#BeginData" << std::endl;

    // Frames
    file << std::hex;
    const auto firstChannel = data.begin();
    for (size_t i = 0; i < firstChannel->second.size(); i++) {
      file << "0x" << std::setw(4) << i;
      for (const auto& [j, channelData] : data) {
        uint16_t sideband = channelData.at(i).valid;
        if (i > 0)
          sideband |= (channelData.at(i).valid and (not channelData.at(i - 1).valid)) << 1;
        if ((i + 1) < channelData.size())
          sideband |= (channelData.at(i).valid and (not channelData.at(i + 1).valid)) << 3;
        file << "    0x" << std::setw(2) << sideband;
        file << " 0x" << std::setw(16) << channelData.at(i).data;
      }
      file << std::endl;
    }
  }

  void writeEMPFile(const BoardData& data, std::ostream& file, const FileFormat format) {
    file << std::setfill('0');

    // Board name/id
    file << "Board CMSSW" << std::endl;

    // Quad/chan header
    file << " Quad/Chan :";
    for (const auto& [i, channelData] : data)
      file << "         q" << std::setw(2) << i / 4 << 'c' << std::setw(1) << i % 4 << "       ";
    file << std::endl;

    // Link header
    file << "      Link :";
    for (const auto& [i, channelData] : data)
      file << "          " << std::setw(3) << i << "        ";
    file << std::endl;

    // Frames
    const auto firstChannel = data.begin();
    for (size_t i = 0; i < firstChannel->second.size(); i++) {
      file << "Frame " << std::setw(4) << i << " :";
      for (const auto& [j, channelData] : data) {
        file << " ";
        //TODO: Add strobe if zero anywhere on channel
        file << "  ";
        file << std::setw(1) << channelData.at(i).valid << "v" << std::setw(16) << std::hex << channelData.at(i).data;
      }
      file << std::endl << std::dec;
    }
  }

  void writeX20File(const BoardData& data, std::ostream& file, const FileFormat format) {
    throw std::runtime_error("Writing X20 file format not yet implemented. Will be done ASAP.");
  }

}  // namespace l1t::demo