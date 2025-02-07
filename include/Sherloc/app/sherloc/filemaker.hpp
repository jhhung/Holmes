#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <future>
#include <spdlog/spdlog.h>
#include <Sherloc/Patient/patient.hpp>
#include <Sherloc/specialcase.hpp>
#include <Sherloc/app/sherloc/sherloc_parameter.hpp>
#include <Sherloc/Attr/utils.hpp>
#include <Sherloc/Attr/inheritance_patterns.hpp>
#include <Sherloc/Attr/rule.hpp>
#include <Sherloc/DB/vep.hpp>
#include <cstdlib>
#include <nlohmann/json.hpp>
#include <spdlog/fmt/ostr.h>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

namespace Sherloc::app::sherloc {

using Path = std::filesystem::path;

class FileMaker
{
public:
  static constexpr auto header_cols = Attr::make_sv_array(
    "chr", "pos", "ref", "alt",
    "HGVSc", "HGVSp", "HGVSg",
    "variant_type", "gene", "transcript", "gene_symbol",
    "rule", "score", "consequence",
    "genotype",
    "clinvar_clnsig", "clinvar_allele_id", "cadd_phred_score", "inheritance_pattern", "dvd_clnsig",
    "extra_info"
  );
  static constexpr auto genotypes = Attr::make_sv_array(
    "unknown", "hete", "homo"
  );
  static constexpr auto consequences = Attr::make_sv_array(
    "pathogenic", "likely_pathogenic", "benign", "likely_benign", "uncertain_significance"
  );
  FileMaker() = default;

  void load_chr(Patient::Patient& patient, std::string_view this_chr){
    using namespace std::literals;
    decltype(auto) para = SherlocParameter::get_paras();
    // TODO: Currently we don't need the INFO part
    auto patient_vcf = DB::HTS_VCF{patient.vcf_file, true, false, false, true};
    auto vcf_status = DB::HTS_VCF::VCF_Status{};
    auto previous_skipped_chr = ""s;
    while((vcf_status = patient_vcf.parse_line()) != DB::HTS_VCF::VCF_Status::VCF_EOF){
      switch (vcf_status) {
        using enum DB::HTS_VCF::VCF_Status;
        case OK:
          break;
        case RECORD_NO_ALT:
          // Assume not to happen
          continue;
        default:
          SPDLOG_ERROR("HTS_VCF parsing error, status code: {}", int(vcf_status));
          exit(1);
      }
      auto& rec = patient_vcf.record;
      try{
        rec.chr = Attr::ChrMap::norm_chr(rec.chr);
      } catch (std::out_of_range& e){
        if(rec.chr != previous_skipped_chr){
          SPDLOG_WARN("Skip seq: '{}', this chromosome is not acceptable", rec.chr);
          previous_skipped_chr = rec.chr;
        }
        continue;
      }

      if(para.detect_sex and rec.chr == "Y"){
        patient.sex = true;
      }

      // TODO: maybe a more elegant way to do this
      if(rec.chr != this_chr){
        continue;
      }

      SherlocMember new_member(rec.chr, rec.pos, rec.ref, rec.alt, rec.genotype);
      SPDLOG_DEBUG("Sher mem GT={}/{}", new_member.genotype[0], new_member.genotype[1]);
      new_member.vcf_format_col = patient_vcf.get_fmt();
      new_member.vcf_id_col = patient_vcf.get_ID();

      // FIXME: there some vcf records have more then one alt
      // e.g. chr3	193643619	.	G	A,T	.	.	.	GT:AD:DP:GQ:PL	1/2:174,123:297:99:2572,0,3905
      // TODO: should make HTS_VCF handle this. And replace the following approach
      if(
        !patient.sher_mems.empty() and 
        patient.sher_mems.back().same_coordinate_as(new_member)
      ) [[unlikely]]
      {
        new_member.next = -1;
        patient.sher_mems.back().next = 1;
      }

      patient.sher_mems.emplace_back(std::move(new_member));
    }
  }

  Path make_damaging_vcf(
    Patient::Patient& patient,
    const Path& vepfile,
    DB::VEP& vep_cache,
    const std::vector<bool>& cache_hit_mask,
    std::string_view this_chr
  ){
    auto damage_vcf_path = vepfile / fmt::format("dmg_for_vep_chr{}_{}.vcf",
      this_chr, patient.name);
    auto ofs = std::ofstream(damage_vcf_path);
    SPDLOG_INFO("[make_damaging_vcf] Saving damage VCF...",
      damage_vcf_path.c_str());
    
    for(auto idx = 0; auto& allele : patient.sher_mems) {
      if(!cache_hit_mask[idx]){ // Those not in cache
        fmt::print(ofs, "chr{}\t{}\t{}\t{}\t{}\t.\t.\n",
          allele.chr, allele.pos, idx, allele.ref, allele.alt);
      }
      ++idx;
    }
    SPDLOG_INFO("[make_damaging_vcf] damage VCF is saved to {}",
      damage_vcf_path.c_str());
    return damage_vcf_path;
  }

  void run_vep(
    Patient::Patient& ze, 
    const Path& vep_output_dir, 
    const DB::VEPRunner& vep_runner,
    const std::vector<std::string>& genes,
    std::string_view this_chr,
    DB::VEP& cache
  ){
    using namespace std::filesystem;
    using namespace DB;

    auto vep_inputfile  = path();
    auto vep_outputfile = vep_output_dir / fmt::format("vep_chr{}_{}.vcf",
      this_chr, ze.name);
    decltype(auto) para = SherlocParameter::get_paras();

    spdlog::stopwatch sw;

    // search variant in cache
    SPDLOG_INFO("Searching cached VEP results ...");
    int total = 0, from_cache = 0;
    sw.reset();
    auto cache_hit_mask = std::vector<bool>(ze.sher_mems.size(), false);
    for(auto& sher_mem : ze.sher_mems){
      if(cache.find(sher_mem).has_value()){
        ++from_cache;
        cache_hit_mask[total] = true;
      }
      ++total;
    }
    SPDLOG_INFO("Take {} sec. chr{}: Cached / Total = {} / {} = {:.2f}%",
      sw, this_chr, from_cache, total, double(from_cache) / (total) * 100.);

    // Task that will parse cache
    auto cache_parsing_task = 
      [&sher_mems = ze.sher_mems, &cache_hit_mask, &cache](bool is_async){
        SPDLOG_INFO("[cache parsing] Parsing cache {}...",
          is_async ? "while vep running" : "");
        auto sw = spdlog::stopwatch{};
        for(auto idx = 0; idx < sher_mems.size(); ++idx){
          if(cache_hit_mask[idx]){
            cache.try_insert_into(sher_mems[idx]);
          }
        }
        SPDLOG_INFO("[cache parsing] Parsing cache takes {} sec.", sw);
      };

    auto all_var_are_in_cache = false;
    // if user want to use existing vep annotation file in `vepfile`, skip vep running
    if(!(para.use_exist_vep_output and std::filesystem::exists(vep_outputfile))){
      vep_inputfile = make_damaging_vcf(ze, vep_output_dir, cache, cache_hit_mask, this_chr);
      if(from_cache < total){ // some variants are not in cache
        SPDLOG_INFO("dmg vcf not empty after cache search, Running VEP to annotate these...");
        auto vep_cmd = vep_runner.make_cmd(
          vep_inputfile, vep_outputfile, genes, para.thread_num, para.grch37);
        SPDLOG_INFO("[run vep] VEP command: {}", vep_cmd);
        SPDLOG_INFO("[run vep] Launch VEP task in background");
        
        sw.reset();
        auto vep_cmd_future = std::async(std::launch::async, std::system, vep_cmd.c_str());

        // parsing cache
        cache_parsing_task(true);

        vep_cmd_future.wait();
        SPDLOG_INFO("VEP runner take {} sec", sw);

        auto vep_cmd_retcode = vep_cmd_future.get();
        if(vep_cmd_retcode != 0){
          SPDLOG_ERROR("[run vep] VEP error, return code = {}", vep_cmd_retcode);
          exit(1);
        }
      }else{
        SPDLOG_INFO("[run vep] Dmg vcf empty after cache search, skip VEP running");
        all_var_are_in_cache = true;
        cache_parsing_task(false);
      }
    }else{
      SPDLOG_INFO("[run vep] Found existing vep output, skipping run_vep.");
      cache_parsing_task(false);
    }
    
    sw.reset();
    if(!all_var_are_in_cache) { // only parse VEP output when the input dmg is not empty
      HTS_VCF vep_output_vcf{vep_outputfile, true, false, true, false};
      VEP::parse_vcf_into(vep_output_vcf, ze.sher_mems);
      SPDLOG_INFO("Parsing VEP output into sherloc member takes {} sec", sw);
    }
    

    // TODO: I think unify all the allele storing convention to VCF style is better
    // consider changing it, otherwise dealing with style conversion is painful
    for(auto& sher_mem : ze.sher_mems){
      if(sher_mem.ref.size() != sher_mem.alt.size()){ // indel vcf style to vep style
        sher_mem.pos++;
        sher_mem.ref = sher_mem.ref.substr(1);
        sher_mem.alt = sher_mem.alt.substr(1);
        if(sher_mem.ref.size() == 0) sher_mem.ref = "-";
        if(sher_mem.alt.size() == 0) sher_mem.alt = "-";
      }
    }
  }

  void run_output( Patient::Patient& ze, std::ofstream& os, const bool output_rule_tag = false) {
    decltype(auto) para = SherlocParameter::get_paras();

    for( auto& sher_mem: ze.sher_mems ) {
      int gt_idx;
      if(sher_mem.genotype[0] == -1 and sher_mem.genotype[1] == -1){
        gt_idx = 0; // unknown genotype
      } else {
        if(sher_mem.genotype[0] != sher_mem.genotype[1])
          gt_idx = 1; // heterozygous genotype
        else
          gt_idx = 2; // homozygous genotype
      }

      // information shared in sherloc_member
      auto sher_mem_info = fmt::format("GnomAD_status={};FORMAT=\"{}\";ID={};{}",
        sher_mem.gnomAD_status,
        sher_mem.vcf_format_col,
        sher_mem.vcf_id_col,
        (sher_mem.gnomAD_status == '2') ? // using 1kg af or not?
          fmt::format("K_AF={};", sher_mem.k_af) :
          ""
      );

      for (int var_idx = 0; var_idx < sher_mem.variants.size(); ++var_idx){
        fmt::print(os,
          "{}\n",
          [&](auto& trans) {
            auto rules = std::vector<Attr::Rule>{};
            std::ranges::set_union(
              sher_mem.group, trans.group,
              std::back_inserter(rules));

            auto tags = std::vector<std::string>{};
            if(output_rule_tag){
              std::ranges::set_union(
                sher_mem.rule_tags, trans.rule_tags,
                std::back_inserter(tags));
            }
            double score = para.filter_rules ?
              para.calculate_score<true>(rules) : 
              para.calculate_score<false>(rules);

            int consequence_idx;
            if      ( score >= 5. )  consequence_idx = 0; // pathogenic
            else if ( score >= 4. )  consequence_idx = 1; // likely pathogenic
            else if ( score <= -5. ) consequence_idx = 2; // benign
            else if ( score <= -4. ) consequence_idx = 3; // likely benign
            else                     consequence_idx = 4; // uncertain significance

            // filtering clinvar & dvd by gene matched or not
            auto clinvar_clnsig = sher_mem.clinvar_valid[var_idx] ?
              sher_mem.clinvar_clnsig :
              "None";
            auto clinvar_allele_id = sher_mem.clinvar_valid[var_idx] ?
              sher_mem.clinvar_allele_id : -1;
            if (sher_mem.clinvar_valid[var_idx]){
              trans.extra_info["Clinvar_geneinfo"] = sher_mem.clinvar_geneinfo;
              trans.extra_info["Clinvar_star"] = std::to_string(sher_mem.clinvar_star);
            }
            auto dvd_clnsig = sher_mem.dvd_valid[var_idx] ?
              sher_mem.dvd_clnsig :
              "None";

            return fmt::format(
              "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}:{}\t{}\t{}{}{}",
              sher_mem.chr, sher_mem.pos, sher_mem.ref, sher_mem.alt, trans.hgvsc, trans.hgvsp, trans.hgvsg,
              fmt::join(trans.type, ";"),
              trans.gene, trans.trans, trans.gene_name,
              fmt::join(rules | std::views::transform(&Attr::Rule::str), ";"),
              score, consequences[consequence_idx],
              genotypes[gt_idx],
              clinvar_clnsig, clinvar_allele_id,
              trans.cadd_phred_score,
              Attr::InheritancePatterns::char2abbreviation(sher_mem.inheritance_patterns[var_idx]), sher_mem.inheritance_pattern_sources[var_idx],
              dvd_clnsig,
              sher_mem_info,
              fmt::join(trans.extra_info | 
                std::views::transform([](const auto& info_pair){
                    return fmt::format("{}={}", info_pair.first, info_pair.second);
                }), ";"),
              output_rule_tag ? fmt::format("\t{}", fmt::join(tags, ";")) : ""
            );
          }(sher_mem.variants[var_idx])
        );
      }
    }
  }
};

}
