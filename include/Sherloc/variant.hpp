#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/serialization/access.hpp>
#include <Sherloc/app/sherloc/sherloc_parameter.hpp>
#include <Sherloc/Attr/utils.hpp>
#include <spdlog/fmt/fmt.h>

namespace Sherloc {

struct Sub_Feature {
  bool has = false;
  bool is_exon = true;
  uint32_t idx = 0;
  uint32_t total = 0;

  Sub_Feature() = default;

  Sub_Feature(bool is_exon, const std::string& idx_total): has(true), is_exon(is_exon) {
    auto split_pos = idx_total.find('/'); // e.g. "37/56"
    idx = std::stoul(idx_total.substr(0, split_pos));
    total = std::stoul(idx_total.substr(split_pos + 1));
  }
};

class Variant {
public:
  

  // vep output STRAND (1/-1)
  //  1 is plus  strand, for feature_strand = true
  // -1 is minus strand, for feature_strand = false
  bool feature_strand = true;

  // plugins
  bool sift = false;
  bool poly = false;
  bool mes_empty = true;
  bool mes = false;
  bool might_escape_nmd = false;
  double mes_score = 0.0;
  double pli = 0.0;
  double revel_score = 0.0;
  double cadd_phred_score = 0.0;

  // sherloc rules
  std::multiset<int> group;
  std::set<std::string> rule_tags;

  // feature
  size_t cds_pos = -1;
  size_t aa_pos = -1;
  std::string gene_name;
  std::string gene;
  std::string trans;
  std::vector<std::string> type;
  std::map<std::string, std::string> extra_info;
  Sub_Feature sub_feat;

  // hgvs
  std::string hgvsc;
  std::string hgvsp;
  std::string hgvsg;

  // amino acid
  std::string codon;

  inline void add_rule(int rule){
    group.emplace(rule);
  }

  inline void add_tag(std::string tag){
    rule_tags.emplace(std::move(tag));
  }

  // After VEP annotation, ensembl id in "Feature" col won't have version number
  // e.g. "ENSTXXXXXXX.V" <- without ".V", but RefSeq ids will have it. e.g. "NM_XXXXXX.V"
  // So we need to remove it to make sure they match DB::GTF keys
  static auto feature_normalize(const std::string& feature){
    if(feature == "-")
      return std::string{""};
    return feature.substr(0, feature.find('.'));
  }

  Variant() = default;

  Variant(
    const Attr::HeaderIndexType& vep_header_index,
    std::string_view vep_record
  ){
    // header: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|REFSEQ_MATCH|REFSEQ_OFFSET|GIVEN_REF|USED_REF|BAM_EDIT|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|HGVSg|CADD_PHRED|CADD_RAW|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|NMD|REVEL|pLI_gene_value
    //         T|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000341065|protein_coding|1/12||ENST00000341065.8:c.4C>T|ENSP00000349216.4:p.His2Tyr|3|4|2|H/Y|Cac/Tac|||1|cds_start_NF|SNV|HGNC|HGNC:28706|||Ensembl||C|C||deleterious_low_confidence(0.02)|benign(0.044)|PANTHER:PTHR12247&PANTHER:PTHR12247:SF67||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000342066|protein_coding|3/14||ENST00000342066.8:c.232C>T|ENSP00000342313.3:p.His78Tyr|322|232|78|H/Y|Cac/Tac|||1||SNV|HGNC|HGNC:28706|||Ensembl||C|C||deleterious(0.04)|possibly_damaging(0.637)|AFDB-ENSP_mappings:AF-Q96NU1-F1.A&PANTHER:PTHR12247&PANTHER:PTHR12247:SF67||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000437963|protein_coding|3/5||ENST00000437963.5:c.232C>T|ENSP00000393181.1:p.His78Tyr|292|232|78|H/Y|Cac/Tac|||1|cds_end_NF|SNV|HGNC|HGNC:28706|||Ensembl||C|C||deleterious(0.04)|possibly_damaging(0.637)|PANTHER:PTHR12247&PANTHER:PTHR12247:SF67&MobiDB_lite:mobidb-lite||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000616016|protein_coding|3/14||ENST00000616016.5:c.769C>T|ENSP00000478421.2:p.His257Tyr|1278|769|257|H/Y|Cac/Tac|||1||SNV|HGNC|HGNC:28706|YES||Ensembl||C|C||deleterious_low_confidence(0.04)|benign(0.332)|PANTHER:PTHR12247&PANTHER:PTHR12247:SF67&MobiDB_lite:mobidb-lite||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000616125|protein_coding|2/11||ENST00000616125.5:c.232C>T|ENSP00000484643.1:p.His78Tyr|232|232|78|H/Y|Cac/Tac|||1|cds_start_NF|SNV|HGNC|HGNC:28706|||Ensembl||C|C||deleterious(0)|possibly_damaging(0.801)|PANTHER:PTHR12247&PANTHER:PTHR12247:SF67||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000617307|protein_coding|2/13||ENST00000617307.5:c.232C>T|ENSP00000482090.2:p.His78Tyr|232|232|78|H/Y|Cac/Tac|||1|cds_start_NF|SNV|HGNC|HGNC:28706|||Ensembl||C|C||deleterious(0.04)|possibly_damaging(0.481)|PANTHER:PTHR12247&PANTHER:PTHR12247:SF67||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000618181|protein_coding|2/10||ENST00000618181.5:c.232C>T|ENSP00000480870.1:p.His78Tyr|232|232|78|H/Y|Cac/Tac|||1|cds_start_NF|SNV|HGNC|HGNC:28706|||Ensembl||C|C||deleterious(0)|possibly_damaging(0.801)|PANTHER:PTHR12247&PANTHER:PTHR12247:SF67&MobiDB_lite:mobidb-lite||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000618323|protein_coding|3/14||ENST00000618323.5:c.769C>T|ENSP00000480678.2:p.His257Tyr|1278|769|257|H/Y|Cac/Tac|||1||SNV|HGNC|HGNC:28706|||Ensembl||C|C||deleterious_low_confidence(0.04)|benign(0.332)|PANTHER:PTHR12247&PANTHER:PTHR12247:SF67&MobiDB_lite:mobidb-lite||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000618779|protein_coding|2/12||ENST00000618779.5:c.232C>T|ENSP00000484256.1:p.His78Tyr|232|232|78|H/Y|Cac/Tac|||1|cds_start_NF|SNV|HGNC|HGNC:28706|||Ensembl||C|C||tolerated(0.05)|possibly_damaging(0.849)|PANTHER:PTHR12247&PANTHER:PTHR12247:SF67||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|missense_variant|MODERATE|SAMD11|ENSG00000187634|Transcript|ENST00000622503|protein_coding|2/13||ENST00000622503.5:c.232C>T|ENSP00000482138.1:p.His78Tyr|232|232|78|H/Y|Cac/Tac|||1|cds_start_NF|SNV|HGNC|HGNC:28706|||Ensembl||C|C||deleterious(0.04)|benign(0.059)|PANTHER:PTHR12247&PANTHER:PTHR12247:SF67||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|missense_variant|MODERATE|SAMD11|148398|Transcript|NM_001385640.1|protein_coding|3/14||NM_001385640.1:c.769C>T|NP_001372569.1:p.His257Tyr|1278|769|257|H/Y|Cac/Tac|||1||SNV|EntrezGene|HGNC:28706|||RefSeq||C|C||deleterious_low_confidence(0.04)|benign(0.332)|||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|missense_variant|MODERATE|SAMD11|148398|Transcript|NM_001385641.1|protein_coding|3/14||NM_001385641.1:c.769C>T|NP_001372570.1:p.His257Tyr|1278|769|257|H/Y|Cac/Tac|||1||SNV|EntrezGene|HGNC:28706|YES||RefSeq||C|C||deleterious_low_confidence(0.04)|benign(0.332)|||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|missense_variant|MODERATE|SAMD11|148398|Transcript|NM_152486.4|protein_coding|3/14||NM_152486.4:c.232C>T|NP_689699.3:p.His78Tyr|322|232|78|H/Y|Cac/Tac|||1||SNV|EntrezGene|HGNC:28706|||RefSeq||C|C||deleterious(0.04)|possibly_damaging(0.637)|||chr1:g.930314C>T|22.4|2.481464|||||0.103|0.00,T|upstream_gene_variant|MODIFIER|LOC107985728|107985728|Transcript|NR_168405.1|lncRNA|||||||||||4710|-1||SNV|EntrezGene||YES||RefSeq||C|C||||||chr1:g.930314C>T|22.4|2.481464||||||
    static constexpr auto pipe_delimiter = Attr::delimiter('|');
    // vep VCF output format will replace ',' with '&'
    static constexpr auto and_delimiter  = Attr::delimiter('&');
    static const std::vector<std::string> extra_info_cols = {
      "SIFT",
      "PolyPhen",
      "REVEL",
      "CADD_PHRED",
      "VARIANT_CLASS",
      "HGNC_ID",
      "MANE_SELECT",
      "MANE_PLUS_CLINICAL",
      "CANONICAL",
      "MaxEntScan_diff",
      "pLI_gene_value"
    };

    std::vector<std::string> entry;
    boost::split(entry, vep_record, pipe_delimiter);
    entry.emplace_back(""); // for unknown col

    // FIXME: if the col index is known, we can avoid using this map approach.
    // this might cause some performance issues
    auto extract_col = [&entry, &vep_header_index](const std::string& key) -> std::string& {
      auto it = vep_header_index.find(key);
      return it != vep_header_index.end() ?
        entry[it->second] :
        entry.back();
    };

    gene = feature_normalize(extract_col("Gene"));
    trans = feature_normalize(extract_col("Feature"));
    boost::split(type, extract_col("Consequence"), and_delimiter);
    auto& aa = extract_col("Amino_acids");
    codon = aa.substr(aa.find('/') + 1);

    // set CDS/AA position (CDS_position	Protein_position)
    {
      auto cds_str = extract_col("CDS_position");
      auto aa_str = extract_col("Protein_position");
      if(!cds_str.empty())
        cds_pos = Attr::parse_vep_pos(cds_str);
      if(!aa_str.empty())
        aa_pos = Attr::parse_vep_pos(aa_str);
    }

    sift = extract_col("SIFT").find("deleterious") != std::string::npos;
    poly = extract_col("PolyPhen").find("probably_damaging") != std::string::npos;

    // REVEL might not have score, so use NAN instead of 0.0
    revel_score = extract_col("REVEL").empty() ? NAN: std::stod(extract_col("REVEL"));

    // CADD might not have score, so use NAN instead of 0.0
    cadd_phred_score = extract_col("CADD_PHRED").empty() ? NAN: std::stod(extract_col("CADD_PHRED"));
    
    mes_empty   = extract_col("MaxEntScan_alt").empty();
    if(!mes_empty){
      double mes_diff   = std::stod(extract_col("MaxEntScan_diff"));
      double mes_ref    = std::stod(extract_col("MaxEntScan_ref"));
      mes       = mes_diff < 0.;
      mes_score = -(mes_diff / mes_ref);
    }
    pli         = extract_col("pLI_gene_value").empty() ? 0. : std::stod(extract_col("pLI_gene_value"));


    feature_strand = (extract_col("STRAND") == "1");
    might_escape_nmd = extract_col("NMD") == "NMD_escaping_variant";

    // set exon / intron positional information
    if(!extract_col("EXON").empty()){
      sub_feat = Sub_Feature{true, extract_col("EXON")};
    } else if(!extract_col("INTRON").empty()) {
      sub_feat = Sub_Feature{false, extract_col("INTRON")};
    }


    // FIXME: these columns are pure strings, is moving them a good approach?
    gene_name   = std::move(extract_col("SYMBOL"));
    hgvsc       = std::move(extract_col("HGVSc"));
    hgvsg       = std::move(extract_col("HGVSg"));
    hgvsp       = std::move(extract_col("HGVSp"));

    // forward these extra information outputs
    for(auto& col : extra_info_cols){
      extra_info.emplace(col, std::move(extract_col(col)));
    }
  }

  static auto make_variants(
    const Attr::HeaderIndexType& vep_header_index,
    std::string_view whole_csq_string
  ){
    std::vector<Variant> variants;
    for(auto&& txp_csq : whole_csq_string | std::views::split(',')){
      variants.emplace_back(
        vep_header_index,std::string_view{std::begin(txp_csq), std::end(txp_csq)});
    }
    return variants;
  }

  void print() {
    fmt::print("variant attrs: {}\t{}\t{}\t{}\t{}\t{}\n",
      gene, trans, hgvsc, hgvsp, hgvsg, codon
    );
    fmt::print("variant rules='{}', types='{}'\n",
      fmt::join(group, ":"),
      fmt::join(type, ",")
    );
  }
};

}
