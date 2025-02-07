#pragma once

#include <vector>
#include <string>
#include <optional>
#include <boost/algorithm/string.hpp>
#include <Sherloc/DB/db.hpp>
#include <Sherloc/DB/vcf.hpp>
#include <Sherloc/Attr/utils.hpp>
#include <Sherloc/sherloc_member.hpp>
#include <Sherloc/Attr/clinical_keywords.hpp>
#include <Sherloc/Attr/inheritance_patterns.hpp>

namespace Sherloc::DB {

class Clinvar
{
public:
  std::string ref = "";
  std::string alt = "";
  std::string clnsig = "None";
  std::string codon = "";

  int64_t allele_id = -1;

  bool onset = false;
  bool severe = false;
  bool consequence = false;
  bool benign = false;

  std::string geneinfo = "";
  int8_t star = 0;
  char adar = 'U';

  HOLMES_SERIALIZE(ar, version){
    ar & ref;
    ar & alt;
    ar & clnsig;
    ar & codon;
    ar & allele_id;
    ar & onset;
    ar & severe;
    ar & consequence;
    ar & benign;
    ar & geneinfo;
    ar & star;
    ar & adar;
  }

  Clinvar() = default;

  Clinvar( const std::string& ref0, const std::string& alt0 )
    :ref(ref0), alt(alt0){}

  Clinvar(HTS_VCF& vcf)
    :ref(vcf.record.ref), alt(vcf.record.alt)
  {
    // ref alt to vep stlye
    std::string clndn, clnvi;
    clndn = vcf.info_str("CLNDN").value_or(std::string{});
    clnsig = vcf.info_str("CLNSIG").value_or(std::string{});
    clnvi = vcf.info_str("CLNVI").value_or(std::string{});
    geneinfo = vcf.info_str("GENEINFO").value_or(std::string{});
    star = status2star(vcf.info_str("CLNREVSTAT").value_or(""));
    std::ranges::for_each(clndn, Attr::as_upper);
    std::ranges::for_each(clnsig, Attr::as_upper);

    for(auto inher_patt : Attr::InheritancePatterns::patterns){
      if(Attr::contains_str(clndn, inher_patt)){
        adar = Attr::InheritancePatterns::pattern2char(inher_patt);
        break;
      }
    }

    onset = Attr::ClinicalKeywords::contains_early_onset_keyword(clndn);
    severe = Attr::ClinicalKeywords::contains_severe_keyword(clndn);
    {
      auto sigs = std::vector<std::string>{};
      boost::split(sigs, clnsig, Attr::delimiter('|'));
      for(auto& sig : sigs){
        if(Attr::ClinicalKeywords::contains_pathogenic_keyword<true>(sig)){
          consequence = true;
          break;
        }
        if(Attr::ClinicalKeywords::contains_benign_keyword<true>(sig)){
          benign = true;
          break;
        }
      }
    }

    std::vector<std::string> variant_identifier;
    boost::split(variant_identifier, clnvi, Attr::delimiter('|'));
    for(auto& vi : variant_identifier){
      if(vi.starts_with("UniProtKB_")){
        auto split_pos = vi.find('#') + 1;
        if(split_pos == 0){ // "#" not found
          split_pos = vi.find(':') + 1;
        }
        codon = vi.substr(split_pos);
        break;
      }
    }

    allele_id = vcf.info_int("ALLELEID").value_or(-1);
  }

  static int8_t status2star(const std::string& status){
    const static auto status_map = std::map<std::string, int8_t>{
      {"no_assertion_criteria_provided", 0},
      {"no_assertion_provided", 0},
      {"no_interpretation_for_the_single_variant", 0},
      {"criteria_provided,_conflicting_interpretations", 1},
      {"criteria_provided,_single_submitter", 1},
      {"criteria_provided,_multiple_submitters,_no_conflicts", 2},
      {"reviewed_by_expert_panel", 3},
      {"practice_guideline", 4},
    };

    auto it = status_map.find(status);
    return (it != status_map.end()) ?
      it->second :
      int8_t{0};
  }

  [[nodiscard]] auto get_genes() const {
    auto genes = std::vector<std::string>{};
    boost::split(genes, geneinfo, Attr::delimiter('|'));
    for(auto& gene : genes){
      // gene = "<genesymbol>:<id>", we want to remove the id
      gene = gene.substr(0, gene.find(':'));
    }
    return genes;
  }
};

class DVD
{
public:
  std::string ref = "";
  std::string alt = "";
  std::string gene_symbol = "";
  std::string clnsig = "None";

  bool onset = false;
  bool severe = false;
  bool consequence = false;
  bool benign = false;

  char adar = 'U';

  HOLMES_SERIALIZE(ar, version){
    ar & ref;
    ar & alt;
    ar & gene_symbol;
    ar & clnsig;
    ar & onset;
    ar & severe;
    ar & consequence;
    ar & benign;
    ar & adar;
  }

  DVD() = default;

  DVD( const std::string& ref0, const std::string& alt0 )
    :ref(ref0), alt(alt0){}

  DVD( HTS_VCF& vcf )
    :ref(vcf.record.ref), alt(vcf.record.alt), gene_symbol(vcf.info_str("GENE").value_or(""))
  {
    auto final_disease = vcf.info_str("FINAL_DISEASE").value_or("");
    auto pathogenicity = vcf.info_str("FINAL_PATHOGENICITY").value_or("");
    std::ranges::for_each(final_disease, [](auto& c){
      Attr::as_upper(c);
      if(c == '_') c = ' ';
    });
    std::ranges::for_each(pathogenicity, Attr::as_upper);
    
    // use dvd info[FINAL_DISEASE] to infer variant AD/AR
    for(auto inher_patt : Attr::InheritancePatterns::patterns){
      if(Attr::contains_str(final_disease, inher_patt)){
        adar = Attr::InheritancePatterns::pattern2char(inher_patt);
        break;
      }
    }

    // use dvd info[FINAL_DISEASE] to infer onset, TODO: check if these keyword really exist
    onset = Attr::ClinicalKeywords::contains_early_onset_keyword(final_disease);

    // use dvd info[FINAL_DISEASE] to infer severe
    severe = Attr::ClinicalKeywords::contains_severe_keyword(final_disease);
    
    consequence = Attr::ClinicalKeywords::contains_pathogenic_keyword<true>(pathogenicity);
    benign = Attr::ClinicalKeywords::contains_benign_keyword<true>(pathogenicity);

    clnsig = pathogenicity;
  }
};

class DataBaseClinvar : public BaseDB {
public:
  // chr -> vec[(pos, Clinvar record)...]
  using PosClinvar = std::pair<size_t, Clinvar>;
  std::map<ChrIndexType, std::vector<PosClinvar>> chr2vec;

  // chr -> {ID of clinvar, index of vec}
  using ID2Index = std::pair<size_t, size_t>;
  std::map<ChrIndexType, std::vector<ID2Index>> clinvar_id2index;


  // Transcript ID -> vec[arr(cds pos, AA pos, ID of clinvar)]
  using PosType = std::array<size_t, 3>;
  std::map<std::string, std::vector<PosType>> txp_map;

  HOLMES_SERIALIZE(ar, _ver) {
    ar & db_version;
    ar & db_build_time;
    ar & chr2vec;
    ar & clinvar_id2index;
    ar & txp_map;
  }

  void from(const Path& filename) override {
    static constexpr auto pipe_delimiter = Attr::delimiter('|');
    // parse Clinvar VCF
    HTS_VCF clinvar_vcf{filename, true, false, true, false};

    db_version = clinvar_vcf
      .get_generic_header_value("fileDate")
      .value_or("None");

    int line_num = 0;
    bool end_vcf = false;

    // parse vep vcf header
    auto desc = clinvar_vcf.get_header_info_description("CSQ");
    SPDLOG_INFO("Clinvar annotated VCF CSQ description: {}", desc);
    auto start_pos = desc.find("Format: ") + 8;
    auto end_pos = desc.find('"', start_pos);
    auto vep_header_index = Attr::make_header_index(
      desc.substr(start_pos, end_pos - start_pos), pipe_delimiter);
    auto index_debug_s = std::string{""};
    for(auto&& [k,v] : vep_header_index){
      index_debug_s += fmt::format("{}:{}\n", k, v);
    }
    SPDLOG_CRITICAL("vep_header_index: \n{}", index_debug_s);

    auto txp_idx        = vep_header_index.at("Feature");
    auto cds_pos_idx    = vep_header_index.at("CDS_position");
    auto aa_idx         = vep_header_index.at("Protein_position");
    auto var_type_idx   = vep_header_index.at("Consequence");
    auto cols = std::vector<std::string>{};

    while (true) {
      switch (clinvar_vcf.parse_line()) {
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
          end_vcf = true;
          break;
        default:
          throw std::runtime_error("unknown status");
      }
      if(end_vcf) break;

      // FIXME: currently DB<Clinvar> store variant allele in VEP style 
      // BUT for indel, position is minus 1 due to some ordering issue
      clinvar_vcf.to_vep_style();
      if(clinvar_vcf.record.ref == "-" or clinvar_vcf.record.alt == "-"){
        --clinvar_vcf.record.pos;
      }

      ChrIndexType chr_idx;
      try{
        chr_idx = Attr::ChrMap::chr2idx(clinvar_vcf.record.chr);
      }catch(std::out_of_range& e){
        SPDLOG_WARN("Skipped unaccepted chr: `{}`", clinvar_vcf.record.chr);
        continue;
      }

      // push ClinVar record
      auto& chromosome = chr2vec[chr_idx];
      chromosome.emplace_back(clinvar_vcf.record.pos, clinvar_vcf);

      // push id2index pair
      auto clinvar_id = std::stoul(clinvar_vcf.get_ID());
      clinvar_id2index[chr_idx].emplace_back(
        clinvar_id, chromosome.size() - 1);

      // parse txp2id pair
      for(auto&& txp_csq : clinvar_vcf
            .info_str("CSQ")
            .value_or("") | std::views::split(','))
      {
        boost::split(
          cols,
          std::string_view{std::begin(txp_csq), std::end(txp_csq)},
          pipe_delimiter);

        // skip non missense variant
        if(cols[var_type_idx].find("missense_variant") == std::string::npos){
          continue;
        }

        auto txp = Variant::feature_normalize(cols[txp_idx]);
        auto& cds_pos_str = cols[cds_pos_idx];
        auto& aa_pos_str = cols[aa_idx];
        if(txp == "" or cds_pos_str == "") continue;

        auto& positions = txp_map[txp];
        size_t cds_pos_num, aa_pos_num;
        try{
          cds_pos_num = Attr::parse_vep_pos(cds_pos_str);
          aa_pos_num = Attr::parse_vep_pos(aa_pos_str);
        }catch(std::exception& e){
          SPDLOG_ERROR("pos, what: {}, cds_pos: {}", e.what(), cds_pos_str);
        }
        positions.emplace_back(std::array{cds_pos_num, aa_pos_num, clinvar_id});
      }

      if(line_num % 500000 == 0){
        SPDLOG_INFO("line {}'s CLNDN: {}, vcf id: {}",
          line_num, clinvar_vcf.info_str("CLNDN").value_or(""), clinvar_id);
      }
      line_num++;
    }

    // check clinvar record positions are sorted
    for(auto& [chr, vec] : chr2vec){
      if(auto til = std::ranges::is_sorted_until(vec, {}, &PosClinvar::first); til != vec.end()){
        auto idx = std::distance(std::begin(vec), til);
        SPDLOG_ERROR("chr{} is sorted til {}", chr, idx);
        SPDLOG_CRITICAL("pos {}, allele id: {}", vec[idx-1].first, vec[idx-1].second.allele_id);
        SPDLOG_CRITICAL("pos {}, allele id: {}", vec[idx].first, vec[idx].second.allele_id);
        SPDLOG_CRITICAL("pos {}, allele id: {}", vec[idx+1].first, vec[idx+1].second.allele_id);
      }
    }

    // sort clinvar id
    for(auto& [chr, id2index]: clinvar_id2index){
      std::ranges::sort(id2index, {}, &ID2Index::first);
    }

    // sort all element in txp_map
    SPDLOG_INFO("Sorting positions...");
    auto cnt = 0;
    for(auto& [txp, positions] : txp_map) {
      std::ranges::sort(positions);
      if(++cnt % 10000 == 0)
        SPDLOG_INFO("{} txps sorted.", cnt);
    }
    SPDLOG_INFO("Done. Total {} txp.", cnt);
  }

  void save(const Path& filename) override {
    this->set_build_time();
    this->log_metadata("DataBaseClinvar");
    save_archive_to(*this, filename);
  }

  void load(const Path& filename) override {
    load_archive_from(*this, filename);
    this->log_metadata("DataBaseClinvar");
  }

  std::optional<Clinvar> find(const std::string& chr0, size_t pos0, const std::string& ref0, const std::string& alt0) {
    auto it_chr = chr2vec.find(Attr::ChrMap::chr2idx(chr0));
    if (it_chr == chr2vec.end())
      return std::nullopt;

    // FIXME: currently Clinvar store variant allele in VEP style 
    // BUT for indel, position is minus 1 due to some ordering issue
    if(ref0 == "-" or alt0 == "-")
      --pos0;

    auto& vec = it_chr->second;

    auto [s_it, e_it] = std::ranges::equal_range(vec, pos0, {}, &PosClinvar::first);
    for(auto it = s_it; it != e_it; ++it){
      auto& clinvar = it->second;
      if(ref0 == clinvar.ref and alt0 == clinvar.alt)
        return clinvar;
    }
    return std::nullopt;
  }

  inline auto find(const SherlocMember& sher_mem) {
    return find(sher_mem.chr, sher_mem.pos, sher_mem.ref, sher_mem.alt);
  }

  Clinvar* find(size_t clinvar_idx){
    for(auto&& [chr, id2index] : clinvar_id2index){
      auto it = std::ranges::lower_bound(id2index, clinvar_idx, {}, &ID2Index::first);
      if (it != id2index.end() and it->first == clinvar_idx){
        return &chr2vec.at(chr)[it->second].second;
      }
    }
    return nullptr;
  }

  static size_t CDS_proj(const PosType& cds_aa_idx){
    return cds_aa_idx.at(0);
  }

  static size_t AA_proj(const PosType& cds_aa_idx){
    return cds_aa_idx.at(1);
  }

  auto find_txp_by(const std::string& txp, std::pair<size_t, size_t> bounds, auto&& proj){
    auto ret = std::vector<Clinvar*>{};
    
    auto txp_it = txp_map.find(txp);
    if(txp_it == txp_map.end()){
      return ret;
    }

    auto clinvar_ids = std::set<size_t>{};
    auto s_it = std::ranges::lower_bound(txp_it->second, bounds.first, {}, proj);
    auto e_it = std::ranges::upper_bound(txp_it->second, bounds.second, {}, proj);
    for(auto it = s_it; it != e_it; ++ it){
      SPDLOG_DEBUG("CDS: {}, AA: {}, CLinvarID: {}",
        it->at(0), it->at(1), it->at(2));
      clinvar_ids.emplace(it->at(2)); // 2 is clinvar id, different txp might repeat so using std::set
    }

    for(auto id : clinvar_ids){
      ret.emplace_back(find(id));
    }
    return ret;
  }
};

class DataBaseDVD : public BaseDB {
public:
  using PosDVD = std::pair<size_t, DVD>;
  std::map<ChrIndexType, std::vector<PosDVD>> db_map;

  HOLMES_SERIALIZE(ar, _ver){
    ar & db_version;
    ar & db_build_time;
    ar & db_map;
  }

  void from(const Path& filename) override {
    HTS_VCF dvd_vcf{filename, true, false, true, false};
    db_version = dvd_vcf.get_generic_header_value("fileDate").value_or("None");

    int line_num = 0;
    HTS_VCF::VCF_Status status;
    while ((status = dvd_vcf.parse_line()) != HTS_VCF::VCF_Status::VCF_EOF) {
      switch (status) {
        using enum HTS_VCF::VCF_Status;
        case OK:
          break;
        case READ_RECORD_FAILED:
          throw std::runtime_error("vcf record unpacking error");
        case UNPACK_FAILED:
          throw std::runtime_error("vcf header parsing error");
        case RECORD_NO_ALT:
          continue;
        default:
          throw std::runtime_error("unknown status");
      }
      // So currently DB<DVD> store variant allele in VEP style
      dvd_vcf.to_vep_style();
      if(dvd_vcf.record.ref == "-" or dvd_vcf.record.alt == "-"){
        --dvd_vcf.record.pos;
      }

      ChrIndexType chr_idx;
      try{
        chr_idx = Attr::ChrMap::chr2idx(dvd_vcf.record.chr);
      }catch(std::out_of_range& e){
        SPDLOG_WARN("Skipped unaccepted chr: `{}`", dvd_vcf.record.chr);
        continue;
      }

      auto& chromosome = db_map[chr_idx];
      chromosome.emplace_back(dvd_vcf.record.pos, dvd_vcf);

      if(line_num % 500000 == 0){
        SPDLOG_INFO("DB<DVD> parsed {} lines.", line_num);
        SPDLOG_INFO("line {}'s FINAL_COMMENTS: {}",
          line_num, dvd_vcf.info_str("FINAL_COMMENTS").value_or("None"));
        dvd_vcf.print();
      }
      line_num++;
    }

    for(auto& [chr, vec] : db_map){
      if(!std::ranges::is_sorted(vec, {}, &PosDVD::first)){
        SPDLOG_ERROR("chr{} is not sorted!", Attr::ChrMap::idx2chr(chr));
        exit(1);
      }
    }
  }

  void save(const Path& filename) override {
    this->set_build_time();
    this->log_metadata("DataBaseDVD");
    save_archive_to(*this, filename);
  }

  void load(const Path& filename) override {
    load_archive_from(*this, filename);
    this->log_metadata("DataBaseDVD");
  }

  std::optional<std::string> get_gene_symbol(const std::string& chr0, size_t pos0) {
    auto it_chr = db_map.find(Attr::ChrMap::chr2idx(chr0));
    if (it_chr == db_map.end())
      return std::nullopt;
    auto&& [_, vec] = *it_chr;
    auto [s_it, e_it] = std::ranges::equal_range(vec, pos0, {}, &PosDVD::first);
    for(auto it = s_it; it != e_it; ++ it){
      if (!it->second.gene_symbol.empty() and it->second.consequence){
        return it->second.gene_symbol;
      }
    }
    return std::nullopt;
  }

  std::optional<DVD> find(const std::string& chr, size_t pos, const std::string& ref, const std::string& alt) {
    auto it_chr = db_map.find(Attr::ChrMap::chr2idx(chr));
    if (it_chr == db_map.end())
      return std::nullopt;

    // FIXME: currently DVD store variant allele in VEP style 
    // BUT for indel, position is minus 1 due to some ordering issue
    if(ref == "-" or alt == "-")
      --pos;

    auto& vec = it_chr->second;
    auto [s_it, e_it] = std::ranges::equal_range(vec, pos, {}, &PosDVD::first);
    for(auto it = s_it; it != e_it; ++ it){
      auto& dvd = it->second;
      if(dvd.ref == ref and dvd.alt == alt){
        return dvd;
      }
    }

    return std::nullopt;
  }

  inline auto find(const SherlocMember& sher_mem) {
    return find(sher_mem.chr, sher_mem.pos, sher_mem.ref, sher_mem.alt);
  }
};

}
