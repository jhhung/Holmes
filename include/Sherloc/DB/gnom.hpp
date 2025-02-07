#pragma once

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <map>
#include <stdexcept>
#include <tuple>
#include <vector>
#include <string>
#include <cstdlib>
#include <Sherloc/DB/exac.hpp>
#include <Sherloc/sherloc_member.hpp>
#include <Sherloc/DB/db.hpp>
#include <Sherloc/DB/vcf.hpp>
#include <omp.h>

namespace Sherloc::DB {

class Gnom_alt
{
  public:
    using PositionType = std::uint32_t;
    using StatusType = char;
    size_t chr;
    size_t pos;
    // for snp
    std::vector< PositionType > db_vec_snp;
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

        ar & db_vec_ins_alt;
        ar & db_vec_del_ref;

        ar & chr;
        ar & pos;
    }

    void reserve(){
        static constexpr size_t snp_capacity = 4000000;
        static constexpr size_t ins_capacity = 300000;
        static constexpr size_t del_capacity = 300000;
        
        db_vec_snp.reserve(snp_capacity);
        db_vec_snp_status.reserve(snp_capacity);

        db_vec_ins.reserve(ins_capacity);
        db_vec_ins_status.reserve(ins_capacity);
        db_vec_ins_alt.reserve(ins_capacity);

        db_vec_del.reserve(del_capacity);
        db_vec_del_status.reserve(del_capacity);
        db_vec_del_ref.reserve(del_capacity);
    }

    void add_allele(HTS_VCF& vcf, bool pass){
        auto pos = (PositionType)vcf.record.pos;
        char status = '0';
        char hom = '0';

        if(pass){
            auto an = vcf.info_int("AN").value();
            auto ac = vcf.info_int("AC").value();
            auto af = an != 0 ? vcf.info_float("AF").value() : 0.f;
            auto nhom = vcf.info_int("nhomalt").value();
            hom += std::min(nhom, (decltype(nhom))2);
            if      ( an < 15000 ) status = '2';
            else if ( af > 0.01  ) status = '7';
            else if ( af > 0.005 ) status = '6';
            else if ( af > 0.003 ) status = '5';
            else if ( af > 0.001 ) status = '4';
            else if ( af >= 0.0005 ) status = '3';
            else status = '8';
        }else{
            status = '1';
        }

        if(vcf.record.ref.size() > vcf.record.alt.size()){ // deletion
            db_vec_del.emplace_back(pos + 1);
            db_vec_del_ref.emplace_back(vcf.record.ref.substr(1));
            db_vec_del_status.emplace_back(compress_status(
                vcf.record.ref[0], status, hom));
        }else if(vcf.record.ref.size() < vcf.record.alt.size()){ // insertion
            db_vec_ins.emplace_back(pos + 1);
            db_vec_ins_alt.emplace_back(vcf.record.alt.substr(1));
            db_vec_ins_status.emplace_back(compress_status(
                vcf.record.alt[0], status, hom));
        }else if(vcf.record.ref.size() == 1){ // snp
            db_vec_snp.emplace_back(pos);
            db_vec_snp_status.emplace_back(compress_status(
                vcf.record.alt[0], status, hom));
        }
    }

    StatusType compress_status(char alt, char status, char hom){
        char compressed_alt;
        switch (alt) {
            case 'A': 
                compressed_alt = 0b00000000; break;
            case 'T': 
                compressed_alt = 0b00000001; break;
            case 'C': 
                compressed_alt = 0b00000010; break;
            case 'G': 
                compressed_alt = 0b00000011; break;
            default:
                throw std::invalid_argument("alt is not goooooood!");
        }
        return compressed_alt 
            | (((status - '0') << 2) & 0b00111100)
            | (((hom    - '0') << 6) & 0b11000000);
    }

    auto decompress_status(StatusType status){
        char alt, st, hom;
        switch (status & 0b00000011) {
            case 0b00000000: 
                alt = 'A'; break;
            case 0b00000001: 
                alt = 'T'; break;
            case 0b00000010: 
                alt = 'C'; break;
            case 0b00000011: 
                alt = 'G'; break;
            default:
                throw std::invalid_argument("status is not goooooood!");
        }
        st  = ((status >> 2) & 0b00001111) + '0';
        hom = ((status >> 6) & 0b00000011) + '0';
        return std::make_tuple(alt, st, hom);
    }

    Exac find_snp( std::uint32_t pos0, char alt0 ){
        auto&& [begin_it, end_it] = std::equal_range(db_vec_snp.begin(), db_vec_snp.end(), pos0);
        Exac exac;
        for(auto it = begin_it; it != end_it; ++it){
            auto idx = std::distance(db_vec_snp.begin(), it);
            std::tie(exac.alt, exac.status, exac.hom) = decompress_status(db_vec_snp_status[idx]);
            if(exac.alt[0] == alt0){
                exac.pos = pos0;
                return exac;
            }
        }
        return {};
    }
    
    Exac find_insertion(std::uint32_t pos0, const std::string& ins_alt){
        auto&& [begin_it, end_it] = std::ranges::equal_range(db_vec_ins, pos0);
        Exac exac;
        char discard;
        for(auto it = begin_it; it != end_it; ++it){
            auto idx = std::distance(std::begin(db_vec_ins), it);
            exac.alt = db_vec_ins_alt[idx];
            if(exac.alt == ins_alt){
                std::tie(discard, exac.status, exac.hom) = decompress_status(db_vec_ins_status[idx]);
                exac.pos = pos0;
                return exac;
            }
        }
        return {};
    }

    Exac find_deletion(std::uint32_t pos0, const std::string& del_ref){
        auto&& [begin_it, end_it] = std::equal_range(db_vec_del.begin(), db_vec_del.end(), pos0);
        Exac exac;
        char discard;
        for(auto it = begin_it; it != end_it; ++it){
            auto idx = std::distance(db_vec_del.begin(), it);
            exac.ref = db_vec_del_ref[idx];
            if(exac.ref == del_ref){
                std::tie(discard, exac.status, exac.hom) = decompress_status(db_vec_del_status[idx]);
                exac.pos = pos0;
                return exac;
            }
        }
        return {};
    }

    inline Exac find( const SherlocMember& sher_mem )
    {
        return find(sher_mem.pos, sher_mem.ref, sher_mem.alt);
    }

    inline Exac find( size_t pos0, const std::string& ref0, const std::string& alt0 )
    {
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

class DataBaseGnomAD : public BaseDB {
public:
  static constexpr size_t chunk_size = 10000000;
  Path gnom_dir;
  std::vector<std::vector<std::string>> db_file;
  Gnom_alt gnom;
  size_t current_chr;
  size_t current_chunk;
  int thread_num = 4;
  bool loaded = false;

  DataBaseGnomAD(const Path& gnom_dir = std::filesystem::temp_directory_path()): gnom_dir(gnom_dir) {}

  inline static auto get_arc_name(size_t current_arc_idx){
    return fmt::format("gnomAD-{}.arc", current_arc_idx);
  }

  static void build_chromosome(const std::string& url, const Path& out_dir) {
    auto gnomad_vcf = HTS_VCF{url, true, true, true, false};
    auto current_arc_idx = size_t{0};
    auto built_gnom = Gnom_alt{};
    built_gnom.reserve();

    auto AC0_filter_idx = gnomad_vcf.filter2id("AC0");
    auto AS_VQSR_filter_idx = gnomad_vcf.filter2id("AS_VQSR");

    auto save_file = [&](){
      auto chr_idx = Attr::ChrMap::chr2idx(gnomad_vcf.record.chr);
      built_gnom.pos = current_arc_idx;
      built_gnom.chr = chr_idx;
      auto chr_dir = out_dir / Attr::ChrMap::idx2chr(chr_idx);
      std::filesystem::create_directories(chr_dir);

      auto output_file = chr_dir / get_arc_name(current_arc_idx);
      save_archive_to(built_gnom, output_file);
      SPDLOG_LOGGER_INFO(spdlog::get("gnomAD-builder"),
        "{} is saved.", output_file.c_str());
    };

    while(true){
      switch (gnomad_vcf.parse_line()) {
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
            save_file();
            return;
        default:
            throw std::runtime_error("unknown status");
      }

      size_t arc_idx = gnomad_vcf.record.pos / chunk_size;
      if(arc_idx != current_arc_idx){
        save_file();
        
        built_gnom = Sherloc::DB::Gnom_alt();
        built_gnom.reserve();
        ++current_arc_idx;
      }

      if(gnomad_vcf.has_filter_id(AC0_filter_idx)) continue;
      built_gnom.add_allele(gnomad_vcf, !gnomad_vcf.has_filter_id(AS_VQSR_filter_idx));
    }
  }

  void from(const Path& urls_filename) override {
    auto is = std::ifstream{urls_filename};
    auto urls = std::vector<std::string>{};
    auto url = std::string{};
    while(is >> url){
      urls.emplace_back(std::move(url));
    }

    auto err_mt_logger = spdlog::stdout_color_mt("gnomAD-builder");

    auto file_size = urls.size();
    omp_set_num_threads(thread_num);
    #pragma omp parallel for
    for(int idx = 0; idx < file_size; ++idx){
      SPDLOG_LOGGER_INFO(err_mt_logger, "Building file: `{}`",
        urls[idx]);
      build_chromosome(urls[idx], gnom_dir);
    }
  }

  void save(const Path& filename) override {
    // pass
  }

  void load(const Path& filename) override {
    gnom_dir = filename;
  }

  Exac find(const std::string& chr0, size_t pos0, const std::string& ref0, const std::string& alt0) {
    auto chr = Attr::ChrMap::chr2idx(chr0);
    auto chr_dir = gnom_dir / Attr::ChrMap::idx2chr(chr);
    if(!std::filesystem::exists(chr_dir)){
      SPDLOG_WARN("chromosome dir: `{}` not exist!", chr_dir.c_str());
      return {};
    }
    size_t arc_idx = pos0 / chunk_size;
    if (chr != current_chr or arc_idx != current_chunk or !loaded) {
      auto arc_file = chr_dir / get_arc_name(arc_idx);
      if (!std::filesystem::exists(arc_file)) {
        return {};
      }
      load_archive_from(gnom, arc_file);
      current_chr = gnom.chr;
      current_chunk = gnom.pos;
      loaded = true;
    }
    return gnom.find(pos0, ref0, alt0);
  }

  inline Exac find(const SherlocMember& sher_mem) {
    return find(sher_mem.chr, sher_mem.pos, sher_mem.ref, sher_mem.alt);
  }
};

}