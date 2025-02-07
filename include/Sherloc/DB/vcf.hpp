#pragma once

#include <vector>
#include <array>
#include <string>
#include <boost/algorithm/string.hpp>
#include <Sherloc/sherloc_member.hpp>
#include <Sherloc/DB/db.hpp>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <spdlog/spdlog.h>
#include <optional>

namespace Sherloc::DB {

class HTS_VCF {
private:
  htsFile *hts_file;
  bcf_hdr_t *vcf_header;
  bcf1_t *vcf_record;

  int unpack_flg = 0;
public:
  enum VCF_Status{
    OK,
    READ_RECORD_FAILED,
    UNPACK_FAILED,
    RECORD_NO_ALT,
    VCF_EOF
  };

  struct VCF_Record{
    std::string chr;
    std::string ref;
    std::string alt;
    int64_t pos = 0;
    std::array<int32_t, 2> genotype;
    std::vector<std::string> other_alt;
    bool phased = false;
  } record;
  
  /**
   * @brief Construct a new hts vcf object
   * 
   * @param vcf_file 
   * @param unpack_alt 
   * @param unpack_info 
   * @param unpack_format 
   */
  HTS_VCF(const Path& vcf_file, bool unpack_alt = true, bool unpack_flt = true, bool unpack_info = true, bool unpack_format = true):
    hts_file(vcf_open(vcf_file.c_str(), "r")), vcf_record(bcf_init())
  {
    SPDLOG_INFO("HTS_VCF: parsing {} header...", vcf_file.c_str());
    if(hts_file == nullptr) {
      throw std::runtime_error("Unable to open file.");
    }
    SPDLOG_DEBUG("HTS_VCF: can open file, try header");
    vcf_header = vcf_hdr_read(hts_file);
    if(vcf_header == nullptr){
      throw std::runtime_error("Unable to read header.");
    }
    SPDLOG_DEBUG("HTS_VCF: can header, done.");

    unpack_flg = 
      (unpack_alt     ? BCF_UN_STR : 0) |
      (unpack_flt     ? BCF_UN_FLT : 0) |
      (unpack_info    ? BCF_UN_INFO: 0) |
      (unpack_format  ? BCF_UN_FMT : 0);
  }

  ~HTS_VCF(){
    bcf_hdr_destroy(vcf_header);
    bcf_destroy(vcf_record); 
    vcf_close(hts_file);
  }

  /**
   * @brief Retrieves the FORMAT field from a VCF record.
   *
   * This function uses the HTSlib library to parse the FORMAT field of a VCF record.
   * It then constructs a string representation of the FORMAT field by joining the field values with colons.
   * If the FORMAT field cannot be parsed, the function returns a string of all emtpy format.
   *
   * @return A string representation of the FORMAT field, order is "GT:AD:DP:GQ:PL"
   */
  inline auto get_fmt() {
    using namespace std::literals;
    static constexpr auto fmt_keys = Attr::make_sv_array("GT","AD","DP","GQ","PL");
    static constexpr auto colon_delimiter = Attr::delimiter(':');
    static const auto empty_ret = "::::"s;
    kstring_t str = KS_INITIALIZE;
    auto vcf_line = ""sv;
    if(vcf_format1(vcf_header, vcf_record, &str) != 0){
      ks_free(&str);
      return empty_ret;
    }else{
      vcf_line = {str.s, str.l - 1}; // ignore '\n'
    }

    for(int i = 0; i < 8; ++i){
      vcf_line.remove_prefix(vcf_line.find('\t') + 1);
    } // assert atleast 8 cols (to INFO) in a line
    
    auto pos = vcf_line.find('\t');
    if(pos == std::string_view::npos){ // no SAMPLE col. FORMAT col is already the last col
      ks_free(&str);
      return empty_ret;
    }

    auto fmt_idx = Attr::make_header_index(vcf_line.substr(0, pos), colon_delimiter);
    vcf_line.remove_prefix(pos + 1);
    auto cols = std::vector<std::string>{};
    boost::split(cols, vcf_line, colon_delimiter);
    
    ks_free(&str);
    return fmt::format("{}",
      fmt::join(fmt_keys |
        std::views::transform([&](auto key){
          auto it = fmt_idx.find(key);
          return it != fmt_idx.end() ?
            cols[it->second] :
            ""s;
        }),
      ":"));
  }

  /**
   * @brief parse a vcf record and return the status of parser 
   * 
   * @return HTS_VCF::VCF_Status 
   */
  auto parse_line(){
    SPDLOG_DEBUG("HTS_VCF: start parse_line");
    if(auto ret = bcf_read(hts_file, vcf_header, vcf_record); ret != 0){
      if(ret == -1)
        return VCF_Status::VCF_EOF;
      else // ret < -1
        return VCF_Status::READ_RECORD_FAILED;
    }
    SPDLOG_DEBUG("HTS_VCF: parse_line try unpack record");
    if(bcf_unpack(vcf_record, unpack_flg) != 0)
      return VCF_Status::UNPACK_FAILED;
    SPDLOG_DEBUG("HTS_VCF: can unpack");
    int *gt = nullptr;
    int n = 0;
    record = {};
    record.chr = std::string{bcf_hdr_id2name(vcf_header, vcf_record->rid)};
    record.pos = vcf_record->pos + 1; // htslib vcf pos is 0 based!!
    record.ref = std::string{vcf_record->d.allele[0]};
    
    // FIXME: doesnot consider about multiple sample vcf
    if(unpack_flg&BCF_UN_STR){ // parse alt allele
      auto n_al = vcf_record->n_allele;
      if(n_al > 2){
        for(int i = 2; i < n_al; ++i){
          record.other_alt.emplace_back(vcf_record->d.allele[i]);
        }
      }else if (n_al < 2){
        return VCF_Status::RECORD_NO_ALT;
      }
      record.alt = std::string{vcf_record->d.allele[1]};
    }
    if(unpack_flg&BCF_UN_FMT){
      bcf_get_genotypes(vcf_header, vcf_record, &gt, &n);
      int pause;
      if(n > 2){ // which is not a normal situation
        std::string fmt = this->get_fmt();
        SPDLOG_CRITICAL("vcf fmt of this record = {}", fmt);
        SPDLOG_CRITICAL("vcf record chr: {}, pos: {}, ref: {}, alt: {}",
          record.chr, record.pos, record.ref, record.alt); 
        std::cin >> pause;
      }
      // assert(n <= 2); // FIXME: assume one sample
      if(n == 0){ // no GT infomation
        record.genotype = {-1, -1};
      }
      else{
        // ex. in htslib parse 1/1 as {4,4}, 1/0 as {4,2}, 1/. as {4,0}
        // TODO: currently we don't take the phased information ('cause we don't have it)
        for(int i = 0; i < n; ++i){
          if(gt[i] == 0)
            record.genotype[i] = -1;
          else{
            record.genotype[i] = (gt[i] >> 1) - 1;
          }
        }
      }
    }
    return VCF_Status::OK;
  }

  /**
   * @brief get integer info by key
   * 
   * @param key 
   * @return std::optional<int64_t>
   */
  inline auto info_int(const char* key){
    int64_t *val = nullptr;
    int n = 0;
    auto status_code = bcf_get_info_int64(vcf_header, vcf_record, key, &val, &n);
    if(status_code < 0){ // failed
      return std::optional<int64_t>(std::nullopt);
    }
    auto ret = val[0];
    if(val) free(val);
    return std::optional<int64_t>(ret);
  }

  /**
   * @brief get floating point info by key
   * 
   * @param key 
   * @return std::optional<float> 
   */
  inline auto info_float(const char* key){
    float *val = nullptr;
    int n = 0;
    auto status_code = bcf_get_info_float(vcf_header, vcf_record, key, &val, &n);
    if(status_code < 0){ // failed
      return std::optional<float>(std::nullopt);
    }
    auto ret = val[0];
    if(val) free(val);
    return std::optional<float>(ret);
  }

  /**
   * @brief get string info by key
   * 
   * @param key 
   * @return std::optional<std::string>
   */
  inline auto info_str(const char* key){
    char *val = nullptr;
    int n = 0;
    auto status_code = bcf_get_info_string(vcf_header, vcf_record, key, &val, &n);
    if(status_code < 0){
      return std::optional<std::string>(std::nullopt);
    }
    auto&& ret = std::string{val};
    if(val) free(val);
    return std::optional<std::string>(std::forward<std::string>(ret));
  }

  /**
   * @brief Convert a VCF filter string to its corresponding ID.
   * @param filter The filter string to convert.
   * @return The ID of the filter, or -1 if not found.
   */
  inline auto filter2id(const char* filter){
    return bcf_hdr_id2int(vcf_header, BCF_DT_ID, filter);;
  }

  /**
   * @brief Check if a VCF record has a given filter ID.
   * @param id The ID of the filter to check.
   * @return True if the record has the filter, false otherwise.
   */
  inline auto has_filter_id(int id){
    for(int i = 0; i < vcf_record->d.n_flt; ++i)
      if(vcf_record->d.flt[i] == id) return true;
    return false;
  }

  /**
   * @brief Get the description of a VCF header INFO field.
   * @param info_name The name of the INFO field.
   * @return The description of the field, or an empty string if not found.
   */
  inline auto get_header_info_description(const char* info_name){
    using namespace std::string_literals;
    auto hrec = bcf_hdr_get_hrec(vcf_header, BCF_HL_INFO, "ID", info_name, "INFO");
    SPDLOG_DEBUG("get hrec, hrec ptr: {:x}", ptrdiff_t(hrec));
    if(!hrec) return ""s;
    int desc_idx = bcf_hrec_find_key(hrec, "Description");
    SPDLOG_DEBUG("get desc_idx, desc_idx: {}", desc_idx);
    if(desc_idx < 0) return ""s;

    return std::string{hrec->vals[desc_idx]};
  }

  /**
   * @brief Get the value of a generic VCF header field.
   * @param header_name The name of the header field.
   * @return The value of the field, or an empty optional if not found.
   */
  inline auto get_generic_header_value(const char* header_name){
    // ##source=IGSRpipeline -> get_header_value("source") should return "IGSRpipeline"
    bcf_hrec_t* hrec = bcf_hdr_get_hrec(vcf_header, BCF_HL_GEN, header_name, nullptr, nullptr);
    if (hrec == nullptr) {
      return std::optional<std::string>{std::nullopt};
    }
    if(hrec->value == nullptr) {
      return std::optional<std::string>{""};
    }
    return std::optional<std::string>{{hrec->value}};
  }

  /**
   * @brief Get the ID of a VCF record.
   * @return The ID of the record.
   */
  inline auto get_ID() {
    return std::string{vcf_record->d.id};
  }

  inline void print() const {
    fmt::print("{}\t{}\t'{}'\t'{}'\tGT={}/{}\n", record.chr, record.pos, record.ref, record.alt, record.genotype[0], record.genotype[1]);
  }

  /**
  * @brief Converts the current VCF record to VEP-style format.
  * 
  * If the variant is indel, adjusts the position and
  * modifies the alleles to match VEP-style format.
  */
  inline void to_vep_style(){
    // TODO: There are some weird variants like this: ref ACGG, alt TACAA
    // How to deal with it?
    if(record.ref.size() != record.alt.size() and 
      (record.ref.size() == 1 or record.alt.size() == 1)){
      ++record.pos;
      if(record.ref.size() == 1){ // insertion
        record.ref[0] = '-';
        record.alt = record.alt.substr(1);
      }else{ // deletion
        record.alt[0] = '-';
        record.ref = record.ref.substr(1);
      }
    }
  }
};

class VCF {
public:
  bool empty = true;

  std::string ref = "";

  std::string alt = "";

  std::array<int32_t, 2> genotype = {-1, -1};

  HOLMES_SERIALIZE(ar, version){
    ar & ref;
    ar & alt;
    ar & genotype;
  }

  VCF(VCF&&) = default;
  VCF& operator =(VCF&&) = default;

  VCF(const VCF&) = default;
  VCF& operator =(const VCF&) = default;

  VCF() = default;

  VCF(
    const std::string& ref0,
    const std::string& alt0,
    const std::array<int32_t, 2>& gt = {-1, -1}) :
      empty(false), ref(ref0), alt(alt0), genotype(gt) {}
};

class DataBaseVcf : public BaseDB {
public:
  std::map<size_t, std::multimap<size_t, VCF>> vcf_map;

  HOLMES_SERIALIZE(ar, _ver) {
    ar & vcf_map;
  }

   void from(const Path& filename) override {
    if(filename.empty()) // if no file provided, skip this
      return;
    auto vcf = HTS_VCF(filename, true, false, false, true);

    while (true) {
      switch (vcf.parse_line()) {
        using enum HTS_VCF::VCF_Status;
        case OK:
          break;
        case READ_RECORD_FAILED:
          throw std::runtime_error("vcf record unpacking error");
        case UNPACK_FAILED:
          throw std::runtime_error("vcf header parsing error");
        case RECORD_NO_ALT:
          continue;
        case VCF_EOF:
          return;
        default:
          throw std::runtime_error("unknown status");
      }
      vcf.to_vep_style();
      auto chr_str = size_t{};
      try{
        chr_str = Attr::ChrMap::chr2idx(vcf.record.chr);
      } catch(std::out_of_range& e) { // skip chr that not in ChrMap
        continue;
      }

      auto& chromosome = vcf_map[chr_str];
      chromosome.emplace(
        vcf.record.pos,
        VCF{vcf.record.ref, vcf.record.alt, vcf.record.genotype});
      
      for(auto& alt : vcf.record.other_alt){
        chromosome.emplace(
          vcf.record.pos,
          VCF{vcf.record.ref, alt, vcf.record.genotype});
      }
    }
  }

  void load(const Path& file_name) override {
    load_archive_from(*this, file_name);
  }
  void save(const Path& file_name) override {
    save_archive_to(*this, file_name);
  }

  VCF* find(const std::string& chr0, size_t pos0, const std::string& ref0, const std::string& alt0) {
    auto chr = Attr::ChrMap::chr2idx(chr0);
    auto it_chr = vcf_map.find(chr);

    if (it_chr == vcf_map.end())
      return nullptr;

    auto [s_it, e_it] = it_chr->second.equal_range(pos0);

    for(auto it = s_it; it != e_it; ++it){
      auto& allele = it->second;
      if(ref0 == allele.ref and alt0 == allele.alt){
        return &allele;
      }
    }
    return nullptr;
  }

  VCF find(const SherlocMember& sher_mem) {
    auto it = find(sher_mem.chr, sher_mem.pos, sher_mem.ref, sher_mem.alt);
    if (it == nullptr)  return {};
    return *it;
  }

  [[nodiscard]] auto empty() const {
    return vcf_map.size() == 0;
  }
};

}
