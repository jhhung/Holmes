#pragma once
#include <Sherloc/sherloc_member.hpp>
#include <Sherloc/specialcase.hpp>
#include <Sherloc/Patient/patient.hpp>
#include <Sherloc/DB/dbset.hpp>
#include <Sherloc/disease_database.hpp>
#include <Sherloc/app/sherloc/sherloc_parameter.hpp>
#include <string>
#include <spdlog/spdlog.h>

namespace Sherloc::app::sherloc {

class Population {
public:
  enum PopFreqTag {
    pf0,
    pf1,
    pf2,
    pf3,
    pf4,
    pf5,
    pf6,
    pf7,
    pf8,
    pf9,
    pf10,
    pf11,
    pf12,
    pf13,
    pf14,
    pf15,
    pf16,
    pf17,
    pf18,
    pf19,
  };
  enum PopHomTag {
    ph0,
    ph1,
    ph2
  };

  void set_hom(
          const SherlocParameter& para
        , SherlocMember& sher_mem
        , DB::DBSet& db, DB::Exac& exac, DB::Coverage& coverage
  ) {
    if (
      coverage.status == '1' or 
      coverage.status == '0' or 
      exac.status <= '2' or // '0', '1', '2'
      exac.hom == '0'
    ) return;

    // gnomAD 1 hom
    if (exac.hom == '1') {
      if (sher_mem.onset and sher_mem.severe){
        sher_mem.add_rule(113);
        sher_mem.add_tag(HOLMES_MAKE_TAG(ph0));
      }
      return;
    }

    // gnomAD 2+ hom
    if (sher_mem.onset and sher_mem.severe) {
      sher_mem.add_rule(116);
      sher_mem.add_tag(HOLMES_MAKE_TAG(ph1));
    } else if (sher_mem.onset or sher_mem.severe) {
      sher_mem.add_rule(117);
      sher_mem.add_tag(HOLMES_MAKE_TAG(ph2));
    }
  }
  
  void set_freq(
          const SherlocParameter& para
        , SherlocMember& sher_mem
        , DB::DBSet& db
        , const Disease& disease, DB::Exac& exac, DB::Coverage& coverage
  ) {
    // exac status 0: absense or has "AC0" filter
    if (exac.status == '0') {
      switch (coverage.status) {
        case '3':
          // coverage status 3: >80% pop has over 20 coverage
          sher_mem.add_rule(135);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf2));
          break;
        case '0':
          // coverage status 0: <10% pop has over 20 coverage
          sher_mem.add_rule(180);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf0));
          break;
        case '1':
        case '2':
          // coverage status 2: pop's average coverage >30
          // coverage status 1: otherwise
          sher_mem.add_rule(179);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf1));
          break;
        default:
          throw std::runtime_error(
            fmt::format("set_freq: exac.status 0, unknown coverage status: {}", coverage.status));
      }
      return;
    }

    // exac status 1: Not pass filter (has "AS_VQSR" filter)
    if (exac.status == '1') {
      sher_mem.add_rule(177);
      sher_mem.add_tag(HOLMES_MAKE_TAG(pf3));
      return;
    }
    
    auto is_AD_or_XLinked = sher_mem.is_autosomal_dominant() or sher_mem.is_x_linked();

    // exac status 2: AN < 15000
    // Use 1kg AF instead
    if (exac.status == '2') {
      auto af = db.db_1kg.find(sher_mem);

      // <0 stands for no result
      if(af < 0.)
        af = 0.;

      sher_mem.k_af = af;

      if (is_AD_or_XLinked) {
        if (af > para.k_ad_high) {
          // Very high (>1.5%)
          sher_mem.add_rule(165);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf12));
        } else if (af > para.k_ad_mid) {
          // High (>0.5%)
          sher_mem.add_rule(166);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf13));
        } else if (af > para.k_ad_low) {
          // Somewhat high (>0.1%)
          sher_mem.add_rule(167);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf14));
        } else {
          // Low
          sher_mem.add_rule(178);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf15));
        }
      } else { // AR / uncertain / complicated
        if (af > para.k_ar_high) {
          // Very high (>3%)
          sher_mem.add_rule(162);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf16));
        } else if (af > para.k_ar_mid) {
          // High (>1%)
          sher_mem.add_rule(163);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf17));
        } else if (af > para.k_ar_low) {
          // Somewhat high (>0.3%)
          sher_mem.add_rule(164);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf18));
        } else {
          // Low
          sher_mem.add_rule(178);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf19));
        }
      }
      return;
    }

    // AN >= 15000
    if (is_AD_or_XLinked) {
      switch (exac.status) {
        case '8':
          // exac status 8: AF < 0.0005 
          // (Actually in sherloc table, it says <8 total alleles. 
          // But we think maybe a ratio is more general. )
          // **Note: In sherloc supplementary figure, they misplaced this rule as EV0161
          sher_mem.add_rule(101);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf7));
          break;
        case '6':
        case '7':
          // exac status 6: 0.01 >= AF > 0.005
          // exac status 7: AF > 0.01
          //  -> AF > 0.005
          sher_mem.add_rule(96);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf4));
          break;
        case '4':
        case '5':
          // exac status 4: 0.003 >= AF > 0.001 
          // exac status 5: 0.005 >= AF > 0.003
          //  -> 0.005 >= AF > 0.001
          sher_mem.add_rule(97);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf5));
          break;
        case '3':
          // exac status 3: 0.001 >= AF >= 0.0005
          sher_mem.add_rule(161);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf6));
          break;
        default:
          throw std::runtime_error(
            fmt::format("set_freq: AN >= 15000 AD/XL, unknown exac status: {}", exac.status));
      }
    } else { // AR / uncertain / complicated
      switch (exac.status) {
        case '8':
          // exac status 8: AF < 0.0005 
          // (Actually in sherloc table, it says <8 total alleles. 
          // But we think maybe a ratio is more general. )
          sher_mem.add_rule(101);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf11));
          break;
        case '7':
          // exac status 7: AF > 0.01
          sher_mem.add_rule(93);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf8));
          break;
        case '6':
        case '5':
          // exac status 6: 0.01 >= AF > 0.005
          // exac status 5: 0.005 >= AF > 0.003
          //  -> 0.01 >= AF > 0.003
          sher_mem.add_rule(94);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf9));
          break;
        case '4':
        case '3':
          // exac status 3: 0.001 >= AF >= 0.0005
          // exac status 4: 0.003 >= AF > 0.001 
          //  -> 0.003 >= AF >= 0.0005
          sher_mem.add_rule(160);
          sher_mem.add_tag(HOLMES_MAKE_TAG(pf10));
          break;
        default:
          throw std::runtime_error(
            fmt::format("set_freq: AN >= 15000 AR, unknown exac status: {}", exac.status));
      }
    }
    return;
  }

  void run(
          Patient::Patient& sher
          , DB::DBSet& db
          , const Disease& disease
          , const std::vector< SpecialCase >& special_cases
  ) {
    decltype(auto) para = SherlocParameter::get_paras();
    for (auto& sher_mem : sher.sher_mems) {
      auto clinvar_opt = db.db_clinvar.find(sher_mem);
      auto dvd_opt = db.db_dvd.find(sher_mem);

      // inheritance patterns for each variants
      auto gene_info_inhe_patts = db.db_gene_info.find(sher_mem);
      SPDLOG_DEBUG("geneinfo inhe: {}", std::string_view{
        std::begin(gene_info_inhe_patts),
        std::end(gene_info_inhe_patts)
      });

      auto dvd_genes = std::set<std::string>{};
      auto dvd_pathogenic = false, dvd_benign = false;
      auto clinvar_genes = std::set<std::string>{};
      auto clinvar_pathogenic = false, clinvar_benign = false;

      if(dvd_opt.has_value()){
        auto& dvd = dvd_opt.value();
        sher_mem.onset = sher_mem.onset or dvd.onset;
        sher_mem.severe = sher_mem.onset or dvd.severe;
        sher_mem.dvd_clnsig = dvd.clnsig;
        if(!dvd.gene_symbol.empty()){
          dvd_genes.emplace(dvd.gene_symbol);
        }

        dvd_pathogenic = dvd.consequence;
        dvd_benign = dvd.benign;
      }

      if(clinvar_opt.has_value()){
        auto& clinvar = clinvar_opt.value();
        sher_mem.onset = sher_mem.onset or clinvar.onset;
        sher_mem.severe = sher_mem.onset or clinvar.severe;
        sher_mem.clinvar_clnsig = clinvar.clnsig;
        sher_mem.clinvar_geneinfo = clinvar.geneinfo;
        sher_mem.clinvar_allele_id = clinvar.allele_id;
        sher_mem.clinvar_star = clinvar.star;
        sher_mem.codon = clinvar.codon;

        auto gene_vec = clinvar.get_genes();
        clinvar_genes = std::set(gene_vec.begin(), gene_vec.end());

        clinvar_pathogenic = clinvar.consequence;
        clinvar_benign = clinvar.benign;
      }

      SPDLOG_DEBUG("DVD GENES: {}", fmt::join(dvd_genes, "|"));
      SPDLOG_DEBUG("CLINVAR GENES: {}", fmt::join(clinvar_genes, "|"));

      // set each variants inheritance patterns and sources
      auto& vars = sher_mem.variants;
      auto& sher_mem_patt = 
        sher_mem.inheritance_patterns = std::vector<char>(vars.size(), 'U');
      auto& sher_mem_src = 
        sher_mem.inheritance_pattern_sources = std::vector<std::string>(vars.size());
      sher_mem.clinvar_valid = std::vector<bool>(vars.size(), false);
      sher_mem.dvd_valid = std::vector<bool>(vars.size(), false);
      for(int idx = 0; idx < vars.size(); ++idx){
        if (dvd_genes.contains(vars[idx].gene_name)){
          sher_mem.dvd_valid[idx] = true;
          if (sher_mem_patt[idx] == 'U') {
            sher_mem_patt[idx] = dvd_opt->adar;
            sher_mem_src[idx] = "DVD";
          }
          if (dvd_pathogenic)
            vars[idx].add_rule(215);
          if (dvd_benign)
            vars[idx].add_rule(216);
        }
        
        if (clinvar_genes.contains(vars[idx].gene_name)) {
          sher_mem.clinvar_valid[idx] = true;
          if (sher_mem_patt[idx] == 'U') {
            sher_mem_patt[idx] = clinvar_opt->adar;
            sher_mem_src[idx] = "ClinVar";
          }
          if (clinvar_pathogenic)
            vars[idx].add_rule(213);
          if (clinvar_benign)
            vars[idx].add_rule(214);
        }
        
        if (gene_info_inhe_patts[idx] != 'U') {
          if (sher_mem_patt[idx] == 'U') {
            sher_mem_patt[idx] = gene_info_inhe_patts[idx];
            sher_mem_src[idx] = "GeneInfo";
          }
        }
        
        // fallback
        if (sher_mem_patt[idx] == 'U') {
          sher_mem_src[idx] = "None";
        }
      }

      // FIXME: use one of the adar for now, but this variable should be deprecated
      // instead, use different adar for different variant
      sher_mem.inheritance_pattern = (
        dvd_opt.has_value() ?
          dvd_opt->adar :
          (clinvar_opt.has_value() ?
            clinvar_opt->adar :
            'U'
          )
      );

      bool is_special_case = false;
      for (const auto & sc : special_cases) {
        if (sher_mem.chr==sc.chr and sher_mem.pos==sc.pos and sher_mem.ref==sc.ref and sher_mem.alt==sc.alt) {
          is_special_case = true;
          break;
        }
      }
      if (is_special_case) continue;
      
      auto exac = db.db_gnom.find(sher_mem);
      auto coverage = db.db_coverage.find(sher_mem);
      sher_mem.gnomAD_status = exac.status;
      set_freq(para, sher_mem, db, disease, exac, coverage);
      set_hom(para, sher_mem, db, exac, coverage);
    }
  }
};

}
