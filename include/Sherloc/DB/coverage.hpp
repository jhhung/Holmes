#pragma once

#include <string>
#include <set>
#include <Sherloc/DB/db.hpp>
#include <Sherloc/DB/hts.hpp>
#include <Sherloc/sherloc_member.hpp>
#include <Sherloc/Attr/utils.hpp>
#include <spdlog/stopwatch.h>

namespace Sherloc::DB {

class Coverage
{
public:
  size_t pos;
  char status = '0';

  Coverage() = default;
  Coverage(const Coverage& cov) = default;
  Coverage(Coverage&& cov) noexcept = default;
  Coverage& operator=(const Coverage& cov) = default;
  Coverage& operator=(Coverage&& cov) noexcept = default;
  
  Coverage(size_t pos0): pos(pos0) {}
  Coverage(size_t pos0, char status0): pos(pos0), status(status0) {}
  Coverage(size_t pos0, double mean, double over_20): pos(pos0) {
    if( over_20 > 0.8 ) status = '3';
    else if( over_20 < 0.1 )  status = '0';
    else if( mean > 30. )  status = '2';
    else status = '1';
  }

  HOLMES_SERIALIZE(ar, version) {
    ar & pos;
    ar & status;
  }

  auto operator<=>(const Coverage& rhs) const {
    return this->pos <=> rhs.pos;
  }
};

class DataBaseCoverage : public BaseDB {
public:
  std::vector<std::set<Coverage>> db_map;

  HOLMES_SERIALIZE(ar, _ver) {
    ar & db_version;
    ar & db_build_time;
    ar & db_map;
  }

  void from(const Path& filename) override {
    static constexpr auto tab_delimiter = Attr::delimiter('\t');
    static constexpr auto version_prefix = std::string_view{"release/"};
    spdlog::stopwatch sw;
    size_t locus_idx, mean_idx, over_20_idx;
    size_t pos_idx = 0;
    bool is_old_format = false;
    auto parse_chr_pos = [&](const std::vector<std::string>& cols, bool is_old_format){
      size_t chr, pos;
      if(is_old_format){
        chr = Attr::ChrMap::chr2idx(cols[locus_idx]);
        pos = std::stoul(cols[pos_idx]);
      }else{
        auto [chr_str, pos_str] = Attr::explode(cols[locus_idx], ':');
        chr = Attr::ChrMap::chr2idx(chr_str);
        pos = std::stoul(pos_str);
      }
      return std::make_pair(chr, pos);
    };

    auto url = std::string{};
    {
      auto is = std::ifstream{filename};
      is >> url;
      auto ver_pos = url.find(version_prefix);
      if(ver_pos != std::string::npos){
        ver_pos += version_prefix.size();
        db_version = url.substr(ver_pos, url.find('/', ver_pos) - ver_pos);
      } else {
        db_version = "None";
      }
    }
    auto hts = HTS_File{url};

    // get header col indices
    int total_cols = 0;
    {
      auto status = hts.parse_line();
      if(status != HTS_File::HTS_Status::OK){
        SPDLOG_ERROR("Can't read `{}`, status code: {}", url, int(status));
        exit(1);
      }
      auto header_index = Attr::make_header_index(hts.line, tab_delimiter);
      is_old_format = !header_index.contains("locus");
      if(is_old_format){ // old format: "chrom pos", new format "locus"
        locus_idx = header_index.at("chrom");
        pos_idx = header_index.at("pos");
      } else {
        locus_idx = header_index.at("locus");
      }
      mean_idx = header_index.at("mean");
      over_20_idx = header_index.at("over_20");
      SPDLOG_INFO("Indices: {}, mean: {}, over20: {} (, pos_idx: {})",
        locus_idx, mean_idx, over_20_idx, pos_idx);
    }

    db_map = std::vector<std::set<Coverage>>(
      Attr::ChrMap::approved_chr.size(), {Coverage{0}});

    size_t pre_chr = 0;
    char pre_status = '0';
    auto cols = std::vector<std::string>{};
    size_t line_num = 0;
    auto status = HTS_File::HTS_Status{};
    while((status = hts.parse_line()) != HTS_File::HTS_EOF){
      if(status == HTS_File::READ_RECORD_FAILED){
        SPDLOG_ERROR("Parsing line num {} failed!", line_num);
        exit(1);
      }
      boost::split(cols, hts.line, tab_delimiter);
      double mean, over_20;
      try{
        mean = std::stod(cols[mean_idx]);
        over_20 = std::stod(cols[over_20_idx]);
      }catch(std::invalid_argument& e){
        SPDLOG_ERROR("Coverage parse double exception: mean = '{}', over_20 = '{}'",
          cols[mean_idx], cols[over_20_idx]);
        exit(1);
      }

      Coverage cov{0, mean, over_20};
      size_t chr;
      try {
        std::tie(chr, cov.pos) = parse_chr_pos(cols, is_old_format);
      } catch (std::invalid_argument& e) {
        SPDLOG_ERROR("line: '{}'", hts.line);
        exit(1);
      }
      if(++line_num % 10000000 == 0){
        SPDLOG_INFO("Parsed {} lines, chr = {}, pos = {}",
          line_num, Attr::ChrMap::idx2chr(chr), cov.pos);
      }
      if(chr != pre_chr){
        SPDLOG_INFO("chr{} done, spend {}s", Attr::ChrMap::idx2chr(pre_chr), sw);
        sw.reset();
        pre_chr = chr;
        pre_status = '0';
      }
      if(cov.status == pre_status)
        continue;
      db_map[chr].emplace(cov);
      pre_status = cov.status;
    }
    SPDLOG_INFO("chr{} done, spend {}s", Attr::ChrMap::idx2chr(pre_chr), sw);
  }

  void save(const Path& filename) override {
    this->set_build_time();
    this->log_metadata("DataBaseCoverage");
    save_archive_to(*this, filename);
  }

  void load(const Path& filename) override {
    load_archive_from(*this, filename);
    this->log_metadata("DataBaseCoverage");
  }

  Coverage find(const std::string& chr0, size_t pos0){
    auto chr_num = Attr::ChrMap::chr2idx(chr0);
    auto it = db_map[chr_num].upper_bound(Coverage{pos0});
    return *std::prev(it);
  }

  Coverage find(const SherlocMember& sher_mem){
    return find(sher_mem.chr, sher_mem.pos);
  }
};

}
