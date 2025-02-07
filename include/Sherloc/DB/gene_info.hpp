#pragma once

#include <fstream>
#include <Sherloc/sherloc_member.hpp>
#include <Sherloc/Attr/inheritance_patterns.hpp>
#include <Sherloc/DB/db.hpp>
#include <boost/algorithm/string.hpp>
#include <set>
#include <spdlog/spdlog.h>
#include <optional>


namespace Sherloc::DB {

class GeneInfo
{
  private:
  public:
    /**
     * @brief All symbols (main or alias) contained in a record of our source file
     * 
     */
    std::set<std::string> symbols = {};

    std::string ensembl_gene_id = "";

    std::string description = "";
    
    /**
     * @brief The inheritance pattern of this gene information record
     * 'D': Autosomal dominant
     * 'R': Autosomal recessive
     * 'X': X-linked
     * 'Y': Y-linked
     * 'U': Unknown
     * TODO: currently we don't consider X-linked dominant / recessive
     */
    char inheritance_pattern = 'U';

    inline void set_inheritance_pattern(std::string_view pattern){
        inheritance_pattern = Attr::InheritancePatterns::pattern2char(pattern);
    }
    
    HOLMES_SERIALIZE(ar, version)
    {
        ar & symbols;
        ar & ensembl_gene_id;
        ar & description;
        ar & inheritance_pattern;
    }

    friend std::ostream& operator<<(std::ostream& os, const GeneInfo& info){
        os << fmt::format("(Symbols: [{}], Ensembl ID: {}, Inheritance: {}, Description: \"{}\")",
            boost::join(info.symbols, ", "), info.ensembl_gene_id, info.inheritance_pattern, info.description);
        return os;
    }
};

class DataBaseGeneInfo : public BaseDB {
  static constexpr auto tab_delimiter = Attr::delimiter('\t');
  static constexpr auto find_inheri_patt = [](const std::string& description){
      auto upper_descr = boost::to_upper_copy(description);
      for(auto inhe_patt : Attr::InheritancePatterns::patterns){
        if(upper_descr.find(inhe_patt) != std::string::npos)
          return inhe_patt;
      }
      return std::string_view{"UNKNOWN"};
    };
public:
  std::map<std::string, std::vector<size_t>> symbol2index;
  std::vector<GeneInfo> gene_info_list;

  HOLMES_SERIALIZE(ar, version) {
    ar & db_version;
    ar & db_build_time;
    ar & symbol2index;
    ar & gene_info_list;
  }

  void parse_omim(const Path& omim_file_name) {
    static constexpr auto version_prefix = std::string_view{"# Generated: "};
    auto omim = std::ifstream{omim_file_name};

    if(omim.is_open()){
      SPDLOG_INFO("Parsing OMIM genemap2 file from {} ...",
        omim_file_name.c_str());
      std::string buffer;
      while(std::getline(omim, buffer)){
        if(buffer.starts_with(version_prefix)){
          db_version += fmt::format("OMIM: [{}], ",
            buffer.substr(version_prefix.size()));
        }
        if(buffer.starts_with("# Chromosome"))
          break;
      }

      // construct tsv header index
      auto header_idx = Attr::make_header_index(buffer.substr(2), tab_delimiter);

      auto entries = std::vector<std::string>{};
      int ad = 0, ar = 0, xl = 0, yl = 0, u = 0, num_records = 0;
      while(std::getline(omim, buffer)){
        if(buffer.starts_with("#")) // reach the end
          break;
        ++num_records;

        // parsing the record
        boost::split(entries, buffer, tab_delimiter);

        GeneInfo info{
          {}, // symbols
          entries[header_idx["Ensembl Gene ID"]], // ensembl gene ID
          entries[header_idx["Phenotypes"]], // description
          'U' // inheritance_pattern
        };

        boost::split(info.symbols, entries[header_idx["Gene Symbols"]],
          boost::is_any_of(", "), boost::token_compress_on);
        
        info.set_inheritance_pattern(find_inheri_patt(info.description));

        switch (info.inheritance_pattern) {
          case 'D': ad++; break;
          case 'R': ar++; break;
          case 'X': xl++; break;
          case 'Y': yl++; break;
          case 'U': u++; break;
          default: break;
        }
        
        gene_info_list.emplace_back(std::move(info));
        auto info_idx = gene_info_list.size() - 1;
        for(auto& symbol : gene_info_list.back().symbols)
          symbol2index[symbol].emplace_back(info_idx);
      }
      SPDLOG_CRITICAL("AD = {}, AR = {}, XLINKED = {}, YLINKED = {}, Unknown = {}",
        ad, ar, xl, yl, u);
      SPDLOG_CRITICAL("Records in OMIM genemap2: {}", num_records);
    }
    else{
      SPDLOG_WARN("OMIM file '{}' is not opened, skip this file", omim_file_name.c_str());
    }
  }

  void parse_ncbi(const Path& ncbi_file_name) {
    auto ncbi = std::ifstream{ncbi_file_name};

    // TODO: set ncbi version
    auto ncbi_ver = std::string{"None"};

    db_version += fmt::format("NCBI Gene: [{}], ", ncbi_ver);

    if(ncbi.is_open()){
      SPDLOG_INFO("Parsing NCBI Gene info file from {} ...",
        ncbi_file_name.c_str());
      std::string buffer;
      ncbi.get(); // ignore the first '#' character
      std::getline(ncbi, buffer);

      // construct tsv header index
      auto header_idx = Attr::make_header_index(buffer, tab_delimiter);

      auto entries = std::vector<std::string>{};
      int ad = 0, ar = 0, xl = 0, yl = 0, u = 0, num_records = 0;
      while(std::getline(ncbi, buffer)){
        ++num_records;
        // parsing the record
        boost::split(entries, buffer, tab_delimiter);

        GeneInfo info{
          {}, // symbols
          "", // ensembl gene ID
          fmt::format("{} {}", // description
            entries[header_idx["description"]], entries[header_idx["Other_designations"]]),
          'U' // inheritance_pattern
        };
        auto& db_refs = entries[header_idx["dbXrefs"]];

        // parsing ensembl gene ID from dbXrefs column
        if(int pos; (pos = db_refs.find("Ensembl:")) != std::string::npos){
          auto last_pos = pos + 7;
          for(; last_pos < db_refs.size(); ++last_pos)
            if(db_refs[last_pos] == '|') break;
          auto ensl_id = db_refs.substr(pos, last_pos - pos);
          if(ensl_id != "-") info.ensembl_gene_id = std::move(ensl_id);
        }

        boost::split(info.symbols, entries[header_idx["Synonyms"]],
          [](auto c){ return c == '|'; });
        
        info.set_inheritance_pattern(find_inheri_patt(info.description));

        switch (info.inheritance_pattern) {
          case 'D': ad++; break;
          case 'R': ar++; break;
          case 'X': xl++; break;
          case 'Y': yl++; break;
          case 'U': u++; break;
          default: break;
        }
        
        gene_info_list.emplace_back(std::move(info));
        auto info_idx = gene_info_list.size() - 1;
        for(auto& symbol : gene_info_list.back().symbols)
          symbol2index[symbol].emplace_back(info_idx);
      }
      SPDLOG_CRITICAL("AD = {}, AR = {}, XLINKED = {}, YLINKED = {}, Unknown = {}",
        ad, ar, xl, yl, u);
      SPDLOG_CRITICAL("Records in NCBI Gene info: {}", num_records);
    }
    else{
      SPDLOG_WARN("NCBI GENE file '{}' is not opened, skip this file",
        ncbi_file_name.c_str());
    }
  }

  void from(const Path& filename) override {
    auto gene_info_input_config = Attr::load_json(filename);
    
    parse_omim(gene_info_input_config["omim"].get<std::string>());
    parse_ncbi(gene_info_input_config["ncbi"].get<std::string>());
  }

  void save(const Path& filename) override {
    this->set_build_time();
    this->log_metadata("DataBaseGeneInfo");
    save_archive_to(*this, filename);
  }

  void load(const Path& filename) override {
    load_archive_from(*this, filename);
    this->log_metadata("DataBaseGeneInfo");
  }


  /**
   * @brief Finds and returns the inheritance patterns of each variant in the given SherlocMember.
   *
   * This function initializes a vector of characters with the same size as the number of variants in the SherlocMember.
   *
   * @param sher_mem A SherlocMember object containing the variants for which to find the inheritance patterns.
   * @return A vector of characters representing the inheritance patterns of each variant in the SherlocMember.
   */
  auto find(const SherlocMember& sher_mem) {
    decltype(auto) vars = sher_mem.variants;
    auto inhe_patts = std::vector<char>(vars.size(), 'U');
    for(int idx = 0; idx < vars.size(); ++idx){
      auto& symbol = vars[idx].gene_name;
      auto available_indices_opt = find(symbol);
      if (!available_indices_opt.has_value()){
        continue;
      }

      auto possible_patts = std::set<char>{};
      for(auto geneinfo_idx : available_indices_opt.value()){
        auto patt = gene_info_list[geneinfo_idx].inheritance_pattern;
        if (patt != 'U'){ // ignore unknown pattern
          possible_patts.emplace(patt);
        }
      }
      
      if (possible_patts.size() == 1){
        inhe_patts[idx] = *possible_patts.begin();
      }
    }
    return inhe_patts;
  }

  std::optional<std::vector<size_t>> find(const std::string& symbol){
    auto it = symbol2index.find((symbol));
    if(it == symbol2index.end()){
      return std::nullopt;
    }
    
    return it->second;
  }
};

}
