#pragma once

#include <Sherloc/DB/uniprot.hpp>
#include <Sherloc/option_parser.hpp>
#include <Sherloc/DB/vep.hpp>
#include <Sherloc/DB/db.hpp>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <spdlog/stopwatch.h>

namespace Sherloc::app::vep_cache_builder {

struct Parameters {
  bool grch37 = false;
  int thread_num;
  std::string vepfile;
  std::string output;
  std::string vepconfig;
};

class GetParameters : public Parameters,
                      public nucleona::app::cli::OptionParser {
public:
  GetParameters(int argc, char const *argv[]) {
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "show help message")
        ("input,i", po::value<std::string>(&vepfile)->
            required(), "input vep annotation file")
        ("output,o", po::value<std::string>(&output)->
            default_value("/tmp/vep_cache/"),
            "output vep cache archive dirname, this dir will end up containing all the `chr.arc` file")
        ("config,c", po::value<std::string>(&vepconfig)->
            default_value(""),
            "vep config file for running vep. If not specified, will assume the `--input` file is already annotated.")
        ("thread_num,t", po::value< int >(&thread_num)->default_value(1), "Thread num (for VEP)")
        ("grch37", po::bool_switch(&grch37), "Use grch37 coordinate")
        ;

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (argc == 1 or vm.count("help")) {
      std::cout << desc << std::endl;
      std::exit(1);
    }

    po::notify(vm);
  }
};

auto run_vep(const GetParameters& args){
  using namespace Sherloc::DB;
  auto runner = VEPRunner(args.vepconfig);
  auto tmp_output_file = Path(fmt::format("{}.vcf", std::tmpnam(nullptr)));
  SPDLOG_INFO("[run vep] temporary output file: {}", tmp_output_file.c_str());

  auto vep_cmd = runner.make_cmd(
    args.vepfile, tmp_output_file, {}, args.thread_num, args.grch37);
  SPDLOG_INFO("[run vep] VEP command: {}", vep_cmd);
  auto ret_code = std::system(vep_cmd.c_str());
  if (ret_code != 0) {
    SPDLOG_ERROR("[run vep] VEP error, return code = {}", ret_code);
    exit(1);
  }
  return tmp_output_file;
}

void Run(const GetParameters &args) {
    // TODO:
    using namespace Sherloc::DB;
    spdlog::stopwatch sw;
    VEP cache;

    Path annotated_vcf{args.vepfile};
    auto to_run_vep = !args.vepconfig.empty();

    if (to_run_vep){
      annotated_vcf = run_vep(args);
    }

    SPDLOG_INFO("Parsing & Saving VEP VCF...");
    sw.reset();
    cache.set_output_dir(args.output);
    cache.from(annotated_vcf);
    SPDLOG_INFO("Parsing & Saving takes {} sec.", sw);

    SPDLOG_INFO("Saving VEP meta data...");
    sw.reset();
    cache.save(args.output);
    SPDLOG_INFO("Saving takes {} sec.", sw);

    // if (to_run_vep){ // clean up intermediate tmpfile
    //   std::filesystem::remove(annotated_vcf);
    // }
    return;
}

} // namespace Sherloc::app::vep_cache_builder
