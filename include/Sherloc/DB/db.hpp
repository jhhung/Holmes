#pragma once

#include <fstream>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/set.hpp>
#include <boost/serialization/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/algorithm/string.hpp>
#include <ctime>
#include <spdlog/spdlog.h>
#include <Sherloc/Attr/clinical_keywords.hpp>
#include <Sherloc/Attr/allele.hpp>

#define HOLMES_SERIALIZE(ar, version) \
  friend class boost::serialization::access; \
  template<class Archive> \
  void serialize(Archive& ar, const unsigned int version)

namespace Sherloc::DB {

using Path = std::filesystem::path;
namespace bios = boost::iostreams;

template <typename T>
concept BoostLoadable = 
  requires(T t, boost::archive::binary_iarchive& arc) {
    { arc & t } -> std::convertible_to<boost::archive::binary_iarchive&>;
  };

template <typename T>
concept BoostSavable = 
  requires(T t, boost::archive::binary_oarchive& arc) {
    { arc & t } -> std::convertible_to<boost::archive::binary_oarchive&>;
  };

template<BoostLoadable DataBase>
inline void load_archive_from(DataBase& db, const std::filesystem::path& file_name){
  // construct boost gzip in stream
  std::ifstream file(file_name, std::ios_base::in | std::ios_base::binary);
  bios::filtering_streambuf<bios::input> fin;
  fin.push(bios::zstd_decompressor());
  fin.push(file);

  // serialize
  auto db_archive = boost::archive::binary_iarchive(fin);
  db_archive & db;
}

template<BoostSavable DataBase>
inline void save_archive_to(DataBase& db, const std::filesystem::path& file_name){
  // construct boost gzip out stream
  std::ofstream file(file_name, std::ios_base::out | std::ios_base::binary);
  bios::filtering_streambuf<bios::output> fout;
  fout.push(bios::zstd_compressor());
  fout.push(file);

  // serialize
  auto db_archive = boost::archive::binary_oarchive(fout);
  db_archive & db;
}

struct BaseDB {
  using TimePointType = long;
  using VersionType = std::string;
  using ChrIndexType = size_t;
  VersionType db_version = "";
  TimePointType db_build_time;

  /**
   * @brief parse / build db from a file/dir
   * 
   */
  virtual void from(const Path&) = 0;

  /**
   * @brief load db from a prebuilt file/dir
   * 
   */
  virtual void load(const Path&) = 0;

  /**
   * @brief save db to a file/dir
   * 
   */
  virtual void save(const Path&) = 0;

  void set_build_time() {
    db_build_time = static_cast<TimePointType>(std::time(nullptr));
  }

  [[nodiscard]] auto get_build_time() const {
    return static_cast<time_t>(db_build_time);
  }
  
  [[nodiscard]] auto get_build_time_str() const {
    static const auto is_blank = boost::is_any_of(" \n\t");
    auto build_time = this->get_build_time();
    auto str = std::string(std::ctime(&build_time));
    while(is_blank(str.back())){
      str.pop_back();
    }
    return str;
  }

  void log_metadata(const std::string& db_class) const {
    SPDLOG_INFO("<{}> Version: [{}], built date: [{}]",
      db_class, this->db_version, this->get_build_time_str());
  }
};

}
