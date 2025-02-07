#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch/catch.hpp>

#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <ranges>
#include <filesystem>
#include <random>

#include <Sherloc/DB/clinvar.hpp>
#include <Sherloc/sherloc_member.hpp>

const auto clinvar_vep_annotated_testdata = std::filesystem::path(DATA_PATH) 
  / "clinvar/clinvar_test.annotated.vcf.gz";

TEST_CASE("Test ClinVar"){
  using namespace Sherloc::DB;
  using namespace std::filesystem;
  auto clinvar_tmp_arc = temp_directory_path() / "holmes_test_clinvar.arc";

  auto clinvar = DataBaseClinvar{};
  clinvar.from(clinvar_vep_annotated_testdata);
  clinvar.save(clinvar_tmp_arc);

  {
    auto clinvar_loaded = DataBaseClinvar{};
    clinvar_loaded.load(clinvar_tmp_arc);
    CHECK(clinvar_loaded.chr2vec.size() == clinvar.chr2vec.size());
    CHECK(clinvar_loaded.clinvar_id2index.size() == clinvar.clinvar_id2index.size());
    CHECK(clinvar_loaded.txp_map.size() == clinvar.txp_map.size());
  }

  for(auto& [chr, vec] : clinvar.chr2vec){
    CHECK(std::ranges::is_sorted(vec, {}, &DataBaseClinvar::PosClinvar::first));
  }

  SECTION("Pathogenic case"){
    // record:
    // 4	87616080	1713002	A	AGCAGCGACAGCAGTGATAGCAGTGACAGCAGTGATAGCAGCGATAGCAGTGACAGCAGCG	.	.	ALLELEID=1770502;CLNDISDB=MONDO:MONDO:0008903,MedGen:C0242379,OMIM:211980|Human_Phenotype_Ontology:HP:0001402,Human_Phenotype_Ontology:HP:0002899,Human_Phenotype_Ontology:HP:0003007,Human_Phenotype_Ontology:HP:0006750,MONDO:MONDO:0007256,MedGen:C2239176,OMIM:114550,Orphanet:88673;CLNDN=Lung_cancer|Hepatocellular_carcinoma;CLNHGVS=NC_000004.12:g.87616080_87616081insGCAGCGACAGCAGTGATAGCAGTGACAGCAGTGATAGCAGCGATAGCAGTGACAGCAGCG;CLNREVSTAT=no_assertion_criteria_provided;CLNSIG=Pathogenic;CLNVC=Insertion;CLNVCSO=SO:0000667;GENEINFO=DSPP:1834;MC=SO:0001820|inframe_indel;ORIGIN=2
    auto patho_case = clinvar.find(
      "4",
      87616081,
      "-",
      "GCAGCGACAGCAGTGATAGCAGTGACAGCAGTGATAGCAGCGATAGCAGTGACAGCAGCG"
    );

    REQUIRE(patho_case.has_value());
    CHECK(patho_case->consequence);
    CHECK(patho_case->allele_id == 1770502);
    CHECK_FALSE(patho_case->benign);
  }

  SECTION("Benign case"){
    // record:
    // 5	143134	2352466	G	A	.	.	ALLELEID=2336495;CLNDISDB=MeSH:D030342,MedGen:C0950123;CLNDN=Inborn_genetic_diseases;CLNHGVS=NC_000005.10:g.143134G>A;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Likely_benign;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;GENEINFO=PLEKHG4B:153478;MC=SO:0001583|missense_variant;ORIGIN=1
    auto benign_case = clinvar.find(
      "5",
      143134,
      "G",
      "A"
    );

    REQUIRE(benign_case.has_value());
    CHECK_FALSE(benign_case->consequence);
    CHECK(benign_case->allele_id == 2336495);
    CHECK(benign_case->benign);
  }

  SECTION("Uncertain significance case"){
    // record:
    // 5	1296371	375479	A	G	.	.	AF_TGP=0.47264;ALLELEID=362290;CLNDISDB=MedGen:C0008707|MedGen:C1840169;CLNDN=Chronic_osteomyelitis|Coronary_artery_disease,_susceptibility_to;CLNHGVS=NC_000005.10:g.1296371A>G;CLNREVSTAT=no_assertion_criteria_provided;CLNSIG=Uncertain_significance|association;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=ClinGen:CA11915807|OMIM:187270.0006;GENEINFO=TERT:7015|LOC110806263:110806263;ORIGIN=513;RS=2735940
    auto uncer_case = clinvar.find(
      "5",
      1296371,
      "A",
      "G"
    );

    REQUIRE(uncer_case.has_value());
    CHECK_FALSE(uncer_case->consequence);
    CHECK(uncer_case->allele_id == 362290);
    CHECK_FALSE(uncer_case->benign);
  }

  SECTION("Not found case"){
    auto not_found_case = clinvar.find(
      "1",
      1000,
      "A",
      "T"
    );

    REQUIRE_FALSE(not_found_case.has_value());
  }

  SECTION("Deletion case"){
    // record:
    // 1	1703608	1810316	CT	C	.	.	ALLELEID=1867317;CLNDISDB=MedGen:CN517202;CLNDN=not_provided;CLNHGVS=NC_000001.11:g.1703609del;CLNREVSTAT=no_assertion_provided;CLNSIG=not_provided;CLNVC=Deletion;CLNVCSO=SO:0000159;GENEINFO=CDK11A:728642;MC=SO:0001589|frameshift_variant,SO:0001619|non-coding_transcript_variant;ORIGIN=1
    auto del_case = clinvar.find(
      "1",
      1703609,
      "T",
      "-"
    );

    REQUIRE(del_case.has_value());
    CHECK(del_case->allele_id == 1867317);
    CHECK_FALSE(del_case->consequence);
    CHECK_FALSE(del_case->benign);
    CHECK(del_case->clnsig == "NOT_PROVIDED");
  }
}