#pragma once

#include <vector>
#include <map>
#include <set>
#include <string>
#include <nlohmann/json.hpp>
#include <Sherloc/Attr/rule.hpp>
#include <Sherloc/Attr/utils.hpp>
#include <filesystem>

namespace Sherloc::app::sherloc {

class SherlocParameter
{
  public:
    int     exac_an         = 15000;
    double  exac_ad_high    = 0.005;
    double  exac_ad_mid     = 0.001;
    int     exac_ac         = 8;
    double  exac_ar_high    = 0.03;
    double  exac_ar_mid     = 0.003;
    double  k_ad_high       = 0.015;
    double  k_ad_mid        = 0.005;
    double  k_ad_low        = 0.001;
    double  k_ar_high       = 0.03;
    double  k_ar_mid        = 0.01;
    double  k_ar_low        = 0.003;
    int     coverage_high   = 8;
    int     coverage_low    = 1;
    int     coverage_in_hom = 3;
    double  yield_rule      = 0.75;

    struct Mutation_Rich_Region{
      int windows = 25;
      int pathogenic_variant_threshold = 4;
    } mut_rich_config;

    // for runtime options
    bool use_exist_vep_output = false;

    // for rule filtering
    bool filter_rules = false;

    // for sex detection
    bool detect_sex = false;

    // assembly version switch
    bool grch37 = false;

    // parallel
    int thread_num = 8;

    std::map< int, double > score_table;

    static SherlocParameter& get_paras()
    {
        static SherlocParameter instance;
        return instance;
    }

    SherlocParameter() {
        score_table[96] = -5;
        score_table[97] = -3;
        score_table[161] = -1;
        score_table[101] = 0.5;
        score_table[135] = 1;
        score_table[93] = -5;
        score_table[94] = -3;
        score_table[160] = -1;
        score_table[165] = -5;
        score_table[166] = -3;
        score_table[167] = -1;
        score_table[178] = 0;
        score_table[162] = -5;
        score_table[163] = -3;
        score_table[164] = -1;
        score_table[177] = 0;
        score_table[179] = 0.5;
        score_table[180] = 0;
        score_table[113] = -1;
        score_table[116] = -4;
        score_table[117] = -2;
        score_table[16] = 5;
        score_table[19] = 2.5;
        score_table[183] = 2;
        score_table[175] = 2.5;
        score_table[194] = 0;
        score_table[17] = 4;
        score_table[196] = 2;
        score_table[184] = 2;
        score_table[198] = 2.5;
        score_table[181] = 2;
        score_table[18] = 4;
        score_table[112] = -3;
        score_table[212] = 4;
        score_table[44] = 1.5;
        score_table[139] = 1.5;
        score_table[115] = 0;
        score_table[172] = 1;
        score_table[185] = 0;
        score_table[114] = 2.5;
        score_table[182] = 2;
        score_table[103] = -2;
        score_table[186] = 0;
        score_table[21] = -5;
        score_table[168] = -2;
        score_table[64] = 5;
        score_table[65] = 5;
        score_table[143] = 3;
        score_table[142] = 5;
        score_table[66] = 3;
        score_table[0] = 2;
        score_table[144] = 2;
        score_table[145] = 2;
        score_table[71] = 2;
        score_table[138] = 4;
        score_table[146] = 2;
        score_table[71] = 2;
        score_table[183] = 2;
        score_table[146] = 2;
        score_table[107] = 0;
        score_table[169] = 1;
        score_table[206] = 2;
        score_table[205] = 4;
        score_table[155] = 1;
        score_table[153] = 2;
        score_table[154] = 1;
        score_table[211] = 0;
        score_table[81] = 1;
        score_table[80] = 2;
        score_table[79] = 3;
        score_table[193] = 1;
        score_table[51] = 1;
        score_table[50] = 2.5;
        score_table[49] = 4;
        score_table[133] = -2.5;
        score_table[60] = -1;
        score_table[132] = -4;
        score_table[61] = -1;
        score_table[84] = -2;
        score_table[130] = -5;
        score_table[134] = -4;
        score_table[53] = -2;
        score_table[140] = -4;
        score_table[141] = -2;
        score_table[129] = -4;
        score_table[33] = -2.5;
        score_table[108] = 0;
        score_table[34] = -1;
        score_table[24] = 1;
        score_table[108] = 0;
        score_table[23] = 2.5;
        score_table[157] = 1;
        score_table[159] = 0;
        score_table[158] = 0.5;
        score_table[26] = 2.5;
        score_table[108] = 0;
        score_table[27] = 1;
        score_table[36] = -2.5;
        score_table[37] = -1;
        score_table[122] = 0.5;
        score_table[109] = 0;
        score_table[126] = -1;
        score_table[121] = 1;
        score_table[127] = -1;
        score_table[125] = -1;
        score_table[119] = 1;
        score_table[118] = -1;
        score_table[187] = 1;
        score_table[188] = 0.5;
        score_table[201] = 0.5;
        score_table[203] = 0;
        score_table[202] = 0.5;
        score_table[191] = -1;
        score_table[213] = 2; // for clinvar pathogenic
        score_table[214] = -2; // for clinvar benign
        score_table[215] = 2; // for dvd pathogenic
        score_table[216] = -2; // for dvd benign
        score_table[217] = 0.5; // for CADD pathogenic
        score_table[218] = -0.5; // for CADD benign
    }

    void load_score_table(const std::filesystem::path& score_table_file){
      auto json_score_table = Attr::load_json(
        score_table_file,
        fmt::format("can't open score table file: {}", score_table_file.c_str())
      );

      for (auto& [rule_id_str, score] : json_score_table.get<std::map<std::string, double>>()){
        auto rule_id = std::stoi(rule_id_str);
        if (!score_table.contains(rule_id)){
          SPDLOG_WARN("score table file contains rule: '{}' that is not presented in default score table, skip it.", rule_id);
          continue;
        }
        auto& default_score = score_table[rule_id];
        if (default_score != score){
          SPDLOG_INFO("default score of rule '{}' is {}, override by provided value: {}",
            rule_id, default_score, score);
        }
        default_score = score;
      }
      SPDLOG_INFO("Done reading score table.");
    }

    template<bool AFScoreFilter = false>
    auto calculate_score(std::vector<Attr::Rule>& rules){
      static std::set<int> rules_need_low_af{
        16, 19, 183, 17, 196, 184, 181, 18, 114, 182, 142, 65, 138, 64
      }; // "... present in the general population at a frequency above the pathogenic range for this gene, do not apply this criteria."

      static std::set<int> rules_exclude_low_af_score{
        19, 183, 17, 196, 184, 181, 18, 114, 182, 142, 65, 138, 64
      };

      static auto check_af_in_pathogenic_range = [](std::vector<Attr::Rule>& rules){
        auto it = std::ranges::find(rules, Attr::Rule(135));
        if(it != rules.end())
          return it;
        return std::ranges::find(rules, Attr::Rule(101));
      };

      auto rule_of_pathogenic_range = check_af_in_pathogenic_range(rules);
      auto af_in_pathogenic_range = (rule_of_pathogenic_range != rules.end());

      double score = 0.;
      bool rules_exclude_low_af_score_applied = false;
      for(auto& r : rules){
        if constexpr (AFScoreFilter){
          // either the rule r doesn't need AF to be in pathogenic range, or need to be in and also be in it
          // then add the score
          if(!rules_need_low_af.contains(r) or af_in_pathogenic_range){
            score += score_table.at(r);
            r.enable();
          }

          // check if any rule r need to exclude low AF rules' score (e.g. EV0135 and EV0101)
          if(rules_exclude_low_af_score.contains(r))
            rules_exclude_low_af_score_applied = true;
        }else{
          score += score_table.at(r);
          r.enable();
        }
      }

      if constexpr (AFScoreFilter){
        if(af_in_pathogenic_range and rules_exclude_low_af_score_applied){
          // remove score of EV0135 or EV0101
          score -= score_table.at(*rule_of_pathogenic_range);
          rule_of_pathogenic_range->disable();
        }
      }
      return score;
    }

    SherlocParameter(const SherlocParameter&) = delete;
    SherlocParameter& operator=(const SherlocParameter&) = delete;
};

}
