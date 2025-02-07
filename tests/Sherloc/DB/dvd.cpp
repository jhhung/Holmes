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

const auto dvd_testdata = std::filesystem::path(DATA_PATH) 
  / "dvd/dvd_test.vcf";

TEST_CASE("Test DVD"){
  using namespace Sherloc::DB;
  using namespace std::filesystem;
  auto dvd_tmp_arc = temp_directory_path() / "holmes_test_dvd.arc";

  auto dvd = DataBaseDVD{};
  dvd.from(dvd_testdata);
  dvd.save(dvd_tmp_arc);

  {
    auto dvd_loaded = DataBaseDVD{};
    dvd_loaded.load(dvd_tmp_arc);
    CHECK(dvd_loaded.db_map.size() == dvd.db_map.size());
    for(auto& [chr, vec] : dvd_loaded.db_map){
      REQUIRE(dvd.db_map.contains(chr));
      CHECK(dvd.db_map.at(chr).size() == vec.size());
    }
  }

  for(auto& [chr, vec] : dvd.db_map){
    CHECK(std::ranges::is_sorted(vec, {}, &DataBaseDVD::PosDVD::first));
  }

  SECTION("Pathogenic case"){
    // record:
    // chr1	6425219	.	G	A	.	PASS	FINAL_COMMENTS=Pathogenicity_is_based_on_ClinVar_submissions._Please_be_aware_that_not_all_submitters_agree_with_this_pathogenicity._;FINAL_DISEASE=.;FINAL_PATHOGENICITY=Pathogenic;FINAL_PMID=.;GENE=ESPN
    auto patho_case = dvd.find(
      "chr1",
      6425219,
      "G",
      "A"
    );

    REQUIRE(patho_case.has_value());
    CHECK(patho_case->consequence);
    CHECK_FALSE(patho_case->benign);
  }

  SECTION("Benign case"){
    // record:
    // chr1	6425205	rs184469259	G	T	.	PASS	FINAL_COMMENTS=This_variant_contains_a_MAF_in_at_least_one_population_that_meets_or_exceeds_our_maximum_cutoff_of_0.005.;FINAL_DISEASE=.;FINAL_PATHOGENICITY=Benign;FINAL_PMID=28492532;GENE=ESPN
    auto benign_case = dvd.find(
      "1",
      6425205,
      "G",
      "T"
    );

    REQUIRE(benign_case.has_value());
    CHECK_FALSE(benign_case->consequence);
    CHECK(benign_case->benign);
  }

  SECTION("Uncertain significance case"){
    // record:
    // chr1	6425207	.	C	G	.	PASS	FINAL_COMMENTS=This_variant_is_a_VUS_because_it_does_not_have_enough_information.;FINAL_DISEASE=.;FINAL_PATHOGENICITY=Unknown_significance;FINAL_PMID=.;GENE=ESPN
    auto uncer_case = dvd.find(
      "1",
      6425207,
      "C",
      "G"
    );

    REQUIRE(uncer_case.has_value());
    CHECK_FALSE(uncer_case->consequence);
    CHECK_FALSE(uncer_case->benign);
  }

  SECTION("Not found case"){
    auto not_found_case = dvd.find(
      "1",
      1000,
      "A",
      "T"
    );

    REQUIRE_FALSE(not_found_case.has_value());
  }

  SECTION("Deletion case"){
    // record:
    // chr1	6426292	rs775045632	AC	A	.	PASS	FINAL_COMMENTS=This_variant_is_a_VUS_because_it_does_not_have_enough_information.;FINAL_DISEASE=.;FINAL_PATHOGENICITY=Unknown_significance;FINAL_PMID=.;GENE=ESPN
    auto del_case = dvd.find(
      "1",
      6426293,
      "C",
      "-"
    );

    REQUIRE(del_case.has_value());
    CHECK_FALSE(del_case->consequence);
    CHECK_FALSE(del_case->benign);
  }

  SECTION("Same position case"){
    // record:
    // chr1	6426487	rs374640448	G	C	.	PASS	FINAL_COMMENTS=This_variant_is_a_VUS_because_it_does_not_have_enough_information.;FINAL_DISEASE=.;FINAL_PATHOGENICITY=Unknown_significance;FINAL_PMID=.;GENE=ESPN
    // chr1	6426487	rs374640448	G	T	.	PASS	FINAL_COMMENTS=This_variant_is_a_VUS_because_it_does_not_have_enough_information.;FINAL_DISEASE=.;FINAL_PATHOGENICITY=Unknown_significance;FINAL_PMID=.;GENE=ESPN
    auto same1 = dvd.find(
      "1",
      6426487,
      "G",
      "C"
    );

    REQUIRE(same1.has_value());
    CHECK_FALSE(same1->consequence);
    CHECK_FALSE(same1->benign);

    auto same2 = dvd.find(
      "1",
      6426487,
      "G",
      "T"
    );

    REQUIRE(same2.has_value());
    CHECK_FALSE(same2->consequence);
    CHECK_FALSE(same2->benign);
  }
}