#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <array>
#include <set>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <Sherloc/app/sherloc/sherloc_parameter.hpp>
#include <Sherloc/variant.hpp>
#include <Sherloc/Attr/allele.hpp>

namespace Sherloc {

class SherlocMember {
public:
  // k genome af
  double k_af = 0.0;

  // Rules used
  // std::vector<int> group;
  std::multiset<int> group;
  std::set<std::string> rule_tags;

  std::vector<Variant> variants;
  std::vector<char> inheritance_patterns;
  std::vector<std::string> inheritance_pattern_sources;

  // chromosome
  std::string chr;

  // position
  size_t pos;

  // reference 
  std::string ref;

  // alteration
  std::string alt;

  // Clinvar significance
  std::string clinvar_clnsig = "None";
  std::string clinvar_geneinfo = "";
  int64_t clinvar_allele_id = -1;
  int8_t clinvar_star = -1;
  std::vector<bool> clinvar_valid;

  // DVD Final pathogenicity
  std::string dvd_clnsig = "None";
  std::vector<bool> dvd_valid;

  int next = 0;

  bool onset = false;

  bool severe = true; // TODO: use user input

  char inheritance_pattern = 'U';

  bool clinical_rule = false;

  std::string codon;

  std::string vcf_format_col = "";
  std::string vcf_id_col = ".";

  // var_type[i][0] = trascript
  // var_type[i][1] = variant type
  // var_type[i][2] = score in variant type condition
  // var_type[i][3+] = rules used
  std::vector< std::vector< std::string > > var_type;

  std::array<int, 2> genotype;

  char gnomAD_status = '0';

  std::string vcf_info = "";

  SherlocMember(const SherlocMember&) = default;
  SherlocMember& operator =(const SherlocMember&) = default;

  SherlocMember(SherlocMember&&) = default;
  SherlocMember& operator =(SherlocMember&&) = default;

  // constructor with subject input vcf
  SherlocMember(const std::string& chr0, size_t pos0, const std::string& ref0, const std::string& alt0,
    std::array<int, 2> genotype0 = {-1, -1}) // default assume 
    : chr{chr0}, pos{pos0}, ref{ref0}, alt{alt0}, genotype(genotype0)
  {}

  inline void add_rule(int rule){
    group.emplace(rule);
  }
  inline void add_tag(std::string tag){
    rule_tags.emplace(std::move(tag));
  }

  auto same_coordinate_as(const SherlocMember& other){
    return 
      other.chr == chr and
      other.pos == pos and
      other.ref == ref and
      other.alt.size() == alt.size();
  }

  void print() {
    fmt::print("{}\t{}\t{}\t{}: rules='{}'\n",
      chr, pos, ref, alt,
      fmt::join(group, ":")
    );
    for (auto& i : variants)
      i.print();
  }

  inline bool adar_unknown(){
    return inheritance_pattern == 'U';
  }

  inline bool is_autosomal_dominant(){
    return inheritance_pattern == 'D';
  }

  inline bool is_autosomal_recessive(){
    return inheritance_pattern == 'R';
  }

  inline bool is_x_linked(){
    return inheritance_pattern == 'X';
  }

  inline bool is_y_linked(){
    return inheritance_pattern == 'Y';
  }

  inline bool gt_unknown(){
    return genotype[0] == -1 or genotype[1] == -1;
  }

  inline bool gt_hetero(){
    return !gt_unknown() and genotype[0] != genotype[1];
  }

  inline bool gt_homo(){
    return !gt_unknown() and genotype[0] == genotype[1];
  }

  inline bool af_above_somewhat_high(){
    static std::set<int> high_af_rules = {
      96, 97, 161, // ad/x high
      93, 94, 160, // ar/xr high
      165, 166, 167, // ad/x (1kg) high
      162, 163, 164 // ar/xr (1kg) high
    };
    return std::ranges::any_of(group,
      [&high_af_rules = high_af_rules](auto rule){
        return high_af_rules.contains(rule);
      });
  }

  [[nodiscard]] inline auto make_id() const {
    return fmt::format("{}_{}_{}_{}",
      Attr::ChrMap::norm_chr(chr), pos, ref, alt);
  }
};



}
