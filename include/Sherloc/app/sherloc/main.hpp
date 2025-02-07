#pragma once

#include <iostream>
#include <fstream>
#include <Sherloc/specialcase.hpp>
#include <Sherloc/option_parser.hpp>
#include <Sherloc/DB/dbset.hpp>
#include <Sherloc/Attr/utils.hpp>
#include <Sherloc/app/sherloc/population.hpp>
#include <Sherloc/app/sherloc/clinical.hpp>
#include <Sherloc/app/sherloc/sherloc_parameter.hpp>
#include <Sherloc/app/sherloc/filemaker.hpp>
#include <Sherloc/app/sherloc/sherloc_consequence.hpp>
#include <Sherloc/app/sherloc/variant_rule.hpp>
#include <Sherloc/app/sherloc/predict.hpp>
#include <Sherloc/Patient/patient.hpp>
#include <Sherloc/Patient/other_patient.hpp>
#include <Sherloc/version.h>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

#define BENCHMARK(op, line, sw) \
  (sw).reset(); line; \
  SPDLOG_INFO("[{}] takes {:.2f} sec.", op, (sw));

namespace Sherloc::app::sherloc {

using Path = std::filesystem::path;

struct Parameters {
  bool disease = false;
  bool consequence = false;
  bool observation = false;
  bool use_exist_vep_output = false;
  bool filter_rules = false;
  bool detect_sex = false;
  bool grch37 = false;
  bool output_rule_tag = false;
  bool version = false;
  bool inspect = false;
  std::string specialcase_file; // txt format, separate each column by '\t'; only chromosome, position, reference, alternative
  std::string json_file;
  std::string score_table_file;
  std::string vepfile;
  std::string vepconfig;
  std::string vepcache;
  std::string dbconfig;
  std::string output;
  std::string db_compression;
  int thread_num;
};

class GetParameters :
  public Parameters,
  public nucleona::app::cli::OptionParser {
public:
  GetParameters(int argc, char const* argv[]) {
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "show help message")
      ("disease,d", po::bool_switch(&disease), "Does disease file need?")
      ("consequence,c", po::bool_switch(&consequence), "Does consequence file need?")
      ("observation,b", po::bool_switch(&observation), "Does observation file need?")
      ("use_exist_vep_output,u", po::bool_switch(&use_exist_vep_output), "Use the existing vep output?")
      ("filter_rules,f", po::bool_switch(&filter_rules), "Apply rule filter")
      ("detect_sex", po::bool_switch(&detect_sex), "Auto detect sex by chrY")
      ("specialcase,s", po::value< std::string >(&specialcase_file)->default_value(""), "special case file")
      ("json,j", po::value< std::string >(&json_file)->default_value("./input.json"), "input json file")
      ("score_table", po::value< std::string >(&score_table_file)->default_value(""), "score table json file that overwrite default table")
      ("vepfile,v", po::value< std::string >(&vepfile)->default_value("./"), "vep intermediate file output path")
      ("vep_config", po::value< std::string >(&vepconfig)->default_value(""),
        "VEP config file, default to ${project dir}/config/vep_config.json")
      ("db_config", po::value< std::string >(&dbconfig)->default_value(""),
        "Holmes database config file, default to ${project dir}/config/db_config.json")
      ("vep_cache", po::value< std::string >(&vepcache)->default_value(""), "VEP cache dir")
      ("thread_num,t", po::value< int >(&thread_num)->default_value(8), "Thread num (mostly for VEP)")
      ("grch37", po::bool_switch(&grch37), "Use grch37 coordinate")
      ("output_rule_tag", po::bool_switch(&output_rule_tag), "Whether to output rule tags for Web App")
      ("output,o", po::value< std::string >(&output)->default_value("result.txt"), "output file")
      ("version,V", po::bool_switch(&version), "Print Holmes Version")
      ("inspect", po::bool_switch(&inspect), "Load DB and print the version & built date. **Note: --vep_config and --db_config must be valid.")
      ;

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (argc == 1 or vm.count("help")) {
      std::cout << desc << std::endl;
      std::exit(1);
    }

    po::notify(vm);
  }
};

void Run(GetParameters& args) {
  using namespace std::literals;
  decltype(auto) para = SherlocParameter::get_paras();
  if(args.version == true){
    fmt::print("Holmes version: {}\n", SHERLOC_VERSION);
    std::exit(0);
  }
  if(args.inspect){
    DB::DBSet::inspect(args.dbconfig);
    return;
  }

  if(!args.score_table_file.empty()){
    para.load_score_table(args.score_table_file);
  }

  // parse vep config file
  auto vep_runner = DB::VEPRunner(
    args.vepconfig.empty() ?
      Path(HOLMES_CONFIG_PATH) / "vep_config.json" :
      Path(args.vepconfig)
  );

  // load database
  auto db = DB::DBSet{
    args.dbconfig.empty() ?
      Path(HOLMES_CONFIG_PATH) / "db_config.json" :
      Path(args.dbconfig),
    vep_runner.get_assembly_file()
  };


  // read input json file
  auto js = Attr::load_json(args.json_file,
    fmt::format("Can't read input json from '{}'", args.json_file));

  auto vep_output_dir = Path{args.vepfile};  

  auto vep_cache    = DB::VEP{};
  auto sher_conseq  = SherlocConsequence{};
  auto op           = Patient::OtherPatient{};
  auto disease      = Disease{};
  auto fm           = FileMaker{};

  if (args.disease == true) {
    disease.empty = false;
  }

  para.use_exist_vep_output = args.use_exist_vep_output;
  para.filter_rules = args.filter_rules;
  para.thread_num = args.thread_num;
  para.detect_sex = args.detect_sex;
  para.grch37 = args.grch37;
  
  std::cout << 
R"(
 _   _       _                     
| | | | ___ | |_ __ ___   ___  ___ 
| |_| |/ _ \| | '_ ` _ \ / _ \/ __|
|  _  | (_) | | | | | | |  __/\__ \
|_| |_|\___/|_|_| |_| |_|\___||___/
)" << std::endl;

  auto sw = spdlog::stopwatch{};

  // load cache
  if(!args.vepcache.empty()){
    BENCHMARK("read vep cache", {
      vep_cache.load(args.vepcache);
    }, sw);
  }else{
    SPDLOG_INFO("No cache loaded.");
  }


  if (args.consequence == true) {
    std::string filename = db.get_base_dir() / "consequence.txt";
    sher_conseq = SherlocConsequence(filename);
  }

  if (args.observation == true) {
    op = Patient::OtherPatient(js["observation"]);
  }

  if (args.disease == true) {
    disease = Disease(js["disease"], db.get_base_dir() / "disease.arc");
  }
  
  SPDLOG_INFO("database read complete.");

  auto population_tree = Population{};
  auto variant_rule_tree = VariantRule{};
  auto clinical_tree = Clinical{};
  auto prediction_tree = Predict{};

  // parse special case file
  auto special_case_list = SpecialCase::load_special_cases(args.specialcase_file);

  auto genes = std::vector<std::string>{};
  {
    std::string gene;
    auto ifs = std::ifstream{js["gene_list_file"]};
    while(ifs >> gene)
      genes.emplace_back(gene);
  }

  // Run pipeline
  for (const auto & patient_json : js["patient"]) {
    auto patient = Patient::Patient{patient_json};

    auto output_file = Path(args.output) / fmt::format("sherloc_{}.txt", patient.name);
    auto os = std::ofstream(output_file);
    if (!os.is_open()) {
      fmt::print(std::cerr, "cannot open file `{}`\n", output_file.c_str());
      exit(1);
    }
    fmt::print(os, "{}{}\n",
      fmt::join(FileMaker::header_cols, "\t"),
      args.output_rule_tag ? "\trule_tag" : ""
    );
    // run by chromosome
    for(auto& this_chr : Attr::ChrMap::approved_chr){
      BENCHMARK(fmt::format("run load chr{}", this_chr),
        fm.load_chr(patient, this_chr), sw);
      if(patient.sher_mems.empty()){
        SPDLOG_WARN("Patient '{}' has no variant on chr{}, skipped",
          patient.name, this_chr);
        continue;
      }
      BENCHMARK("run vep", fm.run_vep(
        patient, vep_output_dir, vep_runner, genes, this_chr, vep_cache), sw);

      BENCHMARK("run population tree", population_tree.run(
        patient, db, disease, special_case_list), sw);
      BENCHMARK("run clinical tree", clinical_tree.run(
        patient, disease, op), sw);
      BENCHMARK("run variant_rule tree", variant_rule_tree.run(
        patient, db, sher_conseq), sw);
      BENCHMARK("run prediction tree", prediction_tree.run(
        patient), sw);
      
      BENCHMARK("write to output file", fm.run_output(
        patient, os, args.output_rule_tag), sw);
      // release processed chr to reduce mem usage
      BENCHMARK(fmt::format("clean up chr{}", this_chr),
        patient.sher_mems.clear(), sw);
    }
    SPDLOG_INFO("Output file path: {}", output_file.c_str());
    SPDLOG_INFO("Patient {} done.", patient.name);
  }

}

}
