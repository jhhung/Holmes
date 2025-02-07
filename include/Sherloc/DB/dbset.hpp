#pragma once

#include <Sherloc/DB/clinvar.hpp>
#include <Sherloc/DB/coverage.hpp>
#include <Sherloc/DB/k.hpp>
#include <Sherloc/DB/gnom.hpp>
#include <Sherloc/DB/fasta.hpp>
#include <Sherloc/DB/gtf.hpp>
#include <Sherloc/DB/uniprot.hpp>
#include <Sherloc/DB/gene_info.hpp>

#define DB_BENCHMARK(db, line, sw) \
  (sw).reset(); line; \
  SPDLOG_INFO("[ DBSet ] Loading <{}> takes {:.2f} sec.", db, (sw));

#define HOLMES_MAKE_TAG(tag) #tag

namespace Sherloc::DB {

struct DBInspector : BaseDB {
  void from(const Path&){}
  void load(const Path&){}
  void save(const Path&){}

  HOLMES_SERIALIZE(ar, _ver){
    ar & db_version;
    ar & db_build_time;
  }

  DBInspector() = delete;
  DBInspector(const Path& filename){
    load_archive_from(*this, filename);
  }
};

class DBSet{
public:
  DataBase1KG       db_1kg;
  DataBaseClinvar   db_clinvar;
  DataBaseDVD       db_dvd;
  DataBaseCoverage  db_coverage;
  DataBaseGnomAD    db_gnom;
  DataBaseGeneInfo  db_gene_info;
  DataBaseFasta     db_fasta;
  DataBaseGTF       db_gtf;
  DataBaseUniprot   db_uniprot;

  Path base_dir;

  static constexpr auto db_names = Attr::make_sv_array(
    "1kg", "clinvar", "dvd", "coverage", "gnom", "gene_info", "gtf", "uniprot"
  );

  static constexpr auto db_paths = Attr::make_sv_array(
    "k_genome.arc", "clinvar.arc", "dvd.arc", "coverage.arc", "gnomAD", "gene_info.arc",
    "homo_gtf.arc", "uniprot.arc"
  );

  static consteval auto get_default(std::string_view db_name){
    auto it = std::ranges::find(db_names, db_name);
    return it != db_names.end() ? 
      db_paths[std::distance(db_names.begin(), it)] : std::string_view{""};
  }

  [[nodiscard]] inline auto get_base_dir() const {
    return base_dir;
  }

  static Path make_base_dir(const Path& base_dir_opt = ""){
    auto ret_base_dir = Path();
    if(!base_dir_opt.empty()){
      ret_base_dir = base_dir_opt;
    }else{
      auto base_dir_default = Path(std::getenv("HOME")) / "holmes_database";
      SPDLOG_WARN("database_dir: is empty. Assume to be `$HOME/holmes_database` = {}",
        base_dir_default.c_str());
      ret_base_dir = base_dir_default;
    }

    Attr::check_exist(ret_base_dir, fmt::format("'{}'", ret_base_dir.c_str()));
    return ret_base_dir;
  }

  static auto make_checked_path(const nlohmann::json& config, const Path& base_dir){
    // check all paths
    auto db_path_config = config["db"]
      .get<std::map<std::string, std::string, std::less<void>>>();
    auto checked_paths = std::map<std::string, Path>{};
    for(auto idx = 0; idx < db_names.size(); ++idx){
      auto db_name = std::string{db_names[idx]};
      auto db_path = 
        // If config has the entry for `db_name` and it's not empty, then we can use it
        // Otherwise, use default relative path from `db_paths`
        Path(db_path_config.contains(db_name) and !db_path_config.at(db_name).empty() ?
          db_path_config.at(db_name) :
          db_paths[idx]
        );

      // If the given path is a relative path, that path should be relative to `base_dir`
      if(db_path.is_relative()){
        db_path = base_dir / db_path;
      }
      Attr::check_exist(db_path,
        fmt::format("db file `{}` not exist!", db_path.c_str()));

      checked_paths.emplace(std::move(db_name), std::move(db_path));
    }
    return checked_paths;
  }

  static void inspect(const Path& config_file){
    auto config = Attr::load_json(config_file,
      fmt::format("DBSet: config_file '{}' can't be opened!", config_file.c_str()));

    auto inpect_base_dir = make_base_dir(config["base_dir"].get<std::string>());

    auto checked_paths = make_checked_path(config, inpect_base_dir);

    for (auto& [name, p] : checked_paths){
      if(std::filesystem::is_directory(p)){ // skip folder DB (like gnomAD)
        continue;
      }

      auto inspector = DBInspector(p);

      // formatting yaml-like version config
      if (name == "coverage"){ // gnomad and coverage are both from gnomAD
        fmt::print("gnomad:\n");
      }else{
        fmt::print("{}:\n",
          name);
      }
      fmt::print("    version: \"{}\"\n", inspector.db_version);
      fmt::print("    build_time: \"{}\"\n", inspector.get_build_time_str());
    }
  }

  DBSet() = default;
  DBSet(const Path& config_file, const std::optional<Path>& assembly_file = std::nullopt){
    auto config = Attr::load_json(config_file,
      fmt::format("DBSet: config_file '{}' can't be opened!", config_file.c_str()));

    base_dir = make_base_dir(config["base_dir"].get<std::string>());

    auto checked_paths = make_checked_path(config, base_dir);

    // load db
    spdlog::stopwatch sw;
    DB_BENCHMARK("1kg", db_1kg.load(checked_paths["1kg"]), sw);
    DB_BENCHMARK("clinvar", db_clinvar.load(checked_paths["clinvar"]), sw);
    DB_BENCHMARK("dvd", db_dvd.load(checked_paths["dvd"]), sw);
    DB_BENCHMARK("coverage", db_coverage.load(checked_paths["coverage"]), sw);
    DB_BENCHMARK("gnom", db_gnom.load(checked_paths["gnom"]), sw);
    DB_BENCHMARK("gene_info", db_gene_info.load(checked_paths["gene_info"]), sw);
    DB_BENCHMARK("fasta", db_fasta.load(
      assembly_file.has_value() ?
        assembly_file.value() :
        checked_paths["fasta"]
    ), sw);
    DB_BENCHMARK("gtf", db_gtf.load(checked_paths["gtf"]), sw);
    DB_BENCHMARK("uniprot", db_uniprot.load(checked_paths["uniprot"]), sw);
  }
};

}
