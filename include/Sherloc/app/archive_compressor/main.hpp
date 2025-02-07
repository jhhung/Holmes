#pragma once

#include <Sherloc/DB/db.hpp>
#include <Sherloc/option_parser.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <fstream>
#include <iostream>
#include <set>
#include <spdlog/stopwatch.h>

namespace Sherloc::app::archive_compressor {

namespace bios = boost::iostreams;

struct Parameters {
  std::string input;
  std::string output;
  std::string zipbackend;
  bool decomp = false;
};

class GetParameters : public Parameters,
                      public nucleona::app::cli::OptionParser {
public:
  GetParameters(int argc, char const *argv[]) {
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "show help message")(
        "input,i", po::value<std::string>(&input)->required(),
        "input arc file")("output,o",
                          po::value<std::string>(&output)->required(),
                          "output arc file")(
        "backend,b", po::value<std::string>(&zipbackend)->default_value("gzip"),
        "compressor backend (gzip, zstd)")("decomp,d", po::bool_switch(&decomp),
                                           "is decompressing");

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (argc == 1 or vm.count("help")) {
      std::cout << desc << std::endl;
      std::exit(1);
    }

    po::notify(vm);
  }
};

void decompress(const DB::Path &input, const DB::Path &output, const std::string& backend) {
  // construct input stream
  auto file = std::ifstream{input, std::ios_base::in | std::ios_base::binary};
  auto fin = bios::filtering_streambuf<bios::input>{};
  if (backend == "gzip") {
    fin.push(bios::gzip_decompressor());
  }else if (backend == "zstd") {
    fin.push(bios::zstd_decompressor());
  }else{
    SPDLOG_ERROR("Unsupported backend {}.", backend);
    exit(1);
  }
  fin.push(file);

  // construct output stream
  auto fout = std::ofstream{output, std::ios_base::out};

  bios::copy(fin, fout);
}

void compress(const DB::Path &input, const DB::Path &output, const std::string& backend) {
  // construct input stream
  auto fin = std::ifstream{input, std::ios_base::in};


  // construct output stream
  auto file = std::ofstream{output, std::ios_base::out | std::ios_base::binary};
  auto fout = bios::filtering_streambuf<bios::output>{};
  if (backend == "gzip") {
    fout.push(bios::gzip_compressor());
  }else if (backend == "zstd") {
    fout.push(bios::zstd_compressor());
  }else{
    SPDLOG_ERROR("Unsupported backend {}.", backend);
    exit(1);
  }
  fout.push(file);

  bios::copy(fin, fout);
}

void Run(GetParameters &&args) {
  using namespace Sherloc::DB;
  spdlog::stopwatch sw;

  if (args.decomp) {
    decompress(args.input, args.output, args.zipbackend);
  }else{
    compress(args.input, args.output, args.zipbackend);
  }

  return;
}

} // namespace Sherloc::app::archive_compressor
