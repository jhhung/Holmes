#pragma once

#include <Sherloc/option_parser.hpp>
#include <Sherloc/DB/vep.hpp>
#include <Sherloc/DB/vcf.hpp>
#include <Sherloc/DB/db.hpp>
#include <spdlog/fmt/fmt.h>
#include <spdlog/fmt/ostr.h>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace Sherloc::app::vep_cache_query {

struct Parameters {
  std::string vep_cache_dir;
  std::string input_vcf;
  std::string output_file;
};

class GetParameters : public Parameters,
                      public nucleona::app::cli::OptionParser {
public:
  GetParameters(int argc, char const *argv[]) {
    namespace po = boost::program_options;
    po::options_description desc("Allowed options");
    desc.add_options()("help,h", "show help message")
        ("input,i", po::value<std::string>(&input_vcf)->
            required(), "input vcf file to query")
        ("output,o", po::value<std::string>(&output_file)->
            required(),
            "output query result tsv file")
        ("vep_cache_dir,v", po::value<std::string>(&vep_cache_dir)->
            required(),
            "vep cache dir")
        ;

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (argc == 1 or vm.count("help")) {
      std::cout << desc << std::endl;
      std::exit(1);
    }

    po::notify(vm);
  }
};

void Run(const GetParameters &args) {
    // TODO:
    using namespace Sherloc::DB;
    using namespace std::literals;
    auto cache = VEP();

    SPDLOG_INFO("Loading VEP cache...");
    cache.load(args.vep_cache_dir);

    auto input_vcf = HTS_VCF(args.input_vcf, true, false, false, false);
    auto vep_header = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|REFSEQ_MATCH|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|HGVSg|CADD_PHRED|CADD_RAW|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|NMD|REVEL|pLI_gene_value"sv;
    auto tsv_header = "chr\tpos\tref\talt\tvep_records";

    auto output_tsv = std::ofstream{args.output_file};
    fmt::print(output_tsv, "# vep header: {}\n", vep_header);
    fmt::print(output_tsv, "{}\n", tsv_header);

    auto vcf_status = HTS_VCF::VCF_Status{};
    auto previous_skipped_chr = ""s;
    int line_id = 0;
    while((vcf_status = input_vcf.parse_line()) != HTS_VCF::VCF_Status::VCF_EOF){
      ++line_id;
      switch (vcf_status) {
        using enum HTS_VCF::VCF_Status;
        case OK:
          break;
        case RECORD_NO_ALT:
          SPDLOG_WARN("vcf record line id: #{} (1-based) has no alt, skip it.", line_id);
          continue;
        default:
          SPDLOG_ERROR("HTS_VCF parsing error, status code: {}", int(vcf_status));
          exit(1);
      }
      auto& rec = input_vcf.record;
      auto normed_chr = ""s;
      try{
        normed_chr = Attr::ChrMap::norm_chr(rec.chr);
      } catch (std::out_of_range& e){
        if(normed_chr != previous_skipped_chr){
          SPDLOG_WARN("Skip seq: '{}', this chromosome is not acceptable", rec.chr);
          previous_skipped_chr = normed_chr;
        }
        continue;
      }


      SherlocMember new_member(normed_chr, rec.pos, rec.ref, rec.alt);
      auto result_it = cache.find(new_member);

      fmt::print(
        output_tsv,
        "{}\t{}\t{}\t{}\t{}\n",
        rec.chr,
        rec.pos,
        rec.ref,
        rec.alt,
        (result_it.has_value() ? *(result_it.value()) : ""s)
      );
    }


    return;
}

} // namespace Sherloc::app::vep_cache_query
