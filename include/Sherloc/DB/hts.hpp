#pragma once

#include <string_view>
#include <filesystem>
#include <htslib/hts.h>

namespace Sherloc::DB {

using Path = std::filesystem::path;

class HTS_File {
protected:
  kstring_t ks = KS_INITIALIZE;
  htsFile *hts_file;
public:
  /**
   * @brief Construct a new hts file reader
   * 
   * @param file 
   */
  HTS_File(const Path& file):
    hts_file(hts_open(file.c_str(), "r"))
  {
    if(hts_file == nullptr) {
      throw std::runtime_error("Unable to open file.");
    }
  }

  enum HTS_Status{
    OK,
    HTS_EOF,
    READ_RECORD_FAILED
  };

  ~HTS_File(){
    ks_free(&ks);
    hts_close(hts_file);
  }

  /**
   * @brief the parsed line
   * 
   * NOTE: This object is just a kstring wrapper, should be carefully accessed
   */
  std::string_view line;

  /**
   * @brief parse an hts line and return the status of parser 
   * 
   * @return HTS_File::HTS_Status 
   */
  auto parse_line(){
    auto status = hts_getline(hts_file, '\n', &ks);
    if(status <= -2){
      return READ_RECORD_FAILED;
    }

    if(status == -1){
      return HTS_EOF;
    }

    line = std::string_view{ks_str(&ks), ks_len(&ks)};
    return OK;
  }
};

}
