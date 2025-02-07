#pragma once

#include <boost/static_assert.hpp>
#include <cstddef>
#include <vector>
#include <string>
#include <boost/algorithm/string.hpp>
#include <Sherloc/sherloc_member.hpp>
#include <cstdlib>
#include <Sherloc/DB/db.hpp>
#include <Sherloc/DB/vcf.hpp>

namespace Sherloc::DB {
class K
{
  public:
    using PositionType = std::uint32_t;
    using StatusType = float;
    size_t chr;

    // for snp
    std::vector< PositionType > db_vec_snp;
    std::vector< char >         db_vec_snp_alt;
    std::vector< StatusType >   db_vec_snp_status;

    // for insertion
    std::vector< PositionType > db_vec_ins;
    std::vector< std::string >  db_vec_ins_alt;
    std::vector< StatusType >   db_vec_ins_status;

    // for deletion
    std::vector< PositionType > db_vec_del;
    std::vector< std::string >  db_vec_del_ref;
    std::vector< StatusType >   db_vec_del_status;
    HOLMES_SERIALIZE(ar, version){
        ar & db_vec_snp_status;
        ar & db_vec_ins_status;
        ar & db_vec_del_status;

        ar & db_vec_snp;
        ar & db_vec_ins;
        ar & db_vec_del;

        ar & db_vec_snp_alt;
        ar & db_vec_ins_alt;
        ar & db_vec_del_ref;

        ar & chr;
    }

    K(size_t chr = 0): chr(chr){ reserve(); }

    void reserve(){
        static constexpr size_t snp_capacity = 7000000;
        static constexpr size_t ins_capacity = 600000;
        static constexpr size_t del_capacity = 600000;
        
        db_vec_snp.reserve(snp_capacity);
        db_vec_snp_status.reserve(snp_capacity);
        db_vec_snp_alt.reserve(snp_capacity);

        db_vec_ins.reserve(ins_capacity);
        db_vec_ins_status.reserve(ins_capacity);
        db_vec_ins_alt.reserve(ins_capacity);

        db_vec_del.reserve(del_capacity);
        db_vec_del_status.reserve(del_capacity);
        db_vec_del_ref.reserve(del_capacity);
    }

    void add_allele(HTS_VCF& vcf){
        auto pos = (PositionType)vcf.record.pos;
        float status = 0.f;

        status = vcf.info_float("AF").value_or(0.f);

        if(vcf.record.ref.size() > vcf.record.alt.size()){ // deletion
            db_vec_del.emplace_back(pos + 1);
            db_vec_del_ref.emplace_back(vcf.record.ref.substr(1));
            db_vec_del_status.emplace_back(status);
        }else if(vcf.record.ref.size() < vcf.record.alt.size()){ // insertion
            db_vec_ins.emplace_back(pos + 1);
            db_vec_ins_alt.emplace_back(vcf.record.alt.substr(1));
            db_vec_ins_status.emplace_back(status);
        }else if(vcf.record.ref.size() == 1){ // snp
            db_vec_snp.emplace_back(pos);
            db_vec_snp_alt.emplace_back(vcf.record.alt[0]);
            db_vec_snp_status.emplace_back(status);
        }else{
            SPDLOG_CRITICAL("WEIRD VARIANT! ref: {}, alt {}", vcf.record.ref, vcf.record.alt);
        }
    }

    float find_snp( std::uint32_t pos0, char alt0 ) const {
        auto&& [begin_it, end_it] = std::ranges::equal_range(db_vec_snp, pos0);
        for(auto it = begin_it; it != end_it; ++it){
            auto idx = std::distance(db_vec_snp.begin(), it);
            if(db_vec_snp_alt[idx] == alt0){
                return db_vec_snp_status[idx];
            }
        }
        return 0.f;
    }
    
    float find_insertion(std::uint32_t pos0, const std::string& ins_alt) const {
        auto&& [begin_it, end_it] = std::ranges::equal_range(db_vec_ins, pos0);
        for(auto it = begin_it; it != end_it; ++it){
            auto idx = std::distance(db_vec_ins.begin(), it);
            if(db_vec_ins_alt[idx] == ins_alt){
                return db_vec_ins_status[idx];
            }
        }
        return 0.f;
    }

    float find_deletion(std::uint32_t pos0, const std::string& del_ref) const {
        auto&& [begin_it, end_it] = std::ranges::equal_range(db_vec_del, pos0);
        for(auto it = begin_it; it != end_it; ++it){
            auto idx = std::distance(db_vec_del.begin(), it);
            if(db_vec_del_ref[idx] == del_ref){
                return db_vec_del_status[idx];
            }
        }
        return 0.f;
    }

    inline float find( const SherlocMember& sher_mem ) const {
        return find(sher_mem.pos, sher_mem.ref, sher_mem.alt);
    }

    inline float find( size_t pos0, const std::string& ref0, const std::string& alt0 )  const {
        if(ref0 != "-" and alt0 != "-"){ // snp
            return find_snp(pos0, alt0[0]);
        }
        else if(ref0 == "-"){ // insertion
            return find_insertion(pos0, alt0);
        }
        // deletion
        return find_deletion(pos0, ref0);
    }

};

class DataBase1KG : public BaseDB {
public:
  std::vector<K> db_map;

  HOLMES_SERIALIZE(ar, _ver) {
    ar & db_version;
    ar & db_build_time;
    ar & db_map;
  }

  DataBase1KG(): db_map(Attr::ChrMap::approved_chr.size()) {}

  void set(const std::string& filename){
    HTS_VCF k_vcf{filename, true, false, true, false};
    int line_num = 0;
    size_t chr = 0;
    db_version = k_vcf
      .get_generic_header_value("fileDate")
      .value_or("None");
    
    auto eof = false;
    while (true) {
      switch (k_vcf.parse_line()) {
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
          eof = true;
          break;
        default:
          throw std::runtime_error("unknown status");
      }
      if(eof){
        break;
      }
      chr = Attr::ChrMap::chr2idx(k_vcf.record.chr);
      db_map[chr].add_allele(k_vcf);
      ++line_num;
      if(line_num % 1000000 == 0){
        SPDLOG_INFO("DB<K> chr{} parsed {} records.",
         k_vcf.record.chr, line_num);
      }
    }
    SPDLOG_INFO("Chr{} Total variant statistics:", k_vcf.record.chr);
    SPDLOG_INFO("SNP count: {}", db_map[chr].db_vec_snp.size());
    SPDLOG_INFO("Insertion count: {}", db_map[chr].db_vec_ins.size());
    SPDLOG_INFO("Deletion count: {}", db_map[chr].db_vec_del.size());
  }

  void from(const Path& filename) override {
    set(filename);
  }

  void from(const std::vector<std::string>& files) {
    for(auto& file : files){
      set(file);
    }
  }

  void load(const Path& filename) override {
    load_archive_from(*this, filename);
    log_metadata("DataBase1KG");
  }

  void save(const Path& filename) override {
    this->set_build_time();
    log_metadata("DataBase1KG");
    save_archive_to(*this, filename);
  }

  inline float find(
    const std::string& chr0,
    size_t pos0,
    const std::string& ref0,
    const std::string& alt0) const {
    size_t chr = Attr::ChrMap::chr2idx(chr0);
    return db_map.at(chr).find(pos0, ref0, alt0);
  }

  inline float find(const SherlocMember& sher_mem) const {
    return find(sher_mem.chr, sher_mem.pos, sher_mem.ref, sher_mem.alt);
  }
};

}
