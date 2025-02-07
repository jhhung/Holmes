#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch/catch.hpp>

#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <ranges>
#include <filesystem>
#include <random>

#include <Sherloc/DB/k.hpp>
#include <Sherloc/sherloc_member.hpp>

auto make_test_k_alt(){
  using namespace Sherloc::DB;
  using namespace std::filesystem;

  auto test_k_file = path(DATA_PATH) / "k_genome/chr1_test.vcf";
  auto k = DataBase1KG{};
  k.from(test_k_file);

  return std::move(k);
}

auto make_multi_chr_k(){
  using namespace Sherloc::DB;
  using namespace std::filesystem;

  auto test_k_file = path(DATA_PATH) / "k_genome/multi_chr_test.vcf";
  auto k = DataBase1KG{};
  k.from(test_k_file);

  return std::move(k);
}

TEST_CASE("Test new K save & load"){
  using namespace Sherloc::DB;
  using namespace std::filesystem;

  auto k_db = make_test_k_alt();
  auto& k = *k_db.db_map.begin();

  auto snp_cnt = k.db_vec_snp.size();
  auto del_cnt = k.db_vec_del.size();
  auto ins_cnt = k.db_vec_ins.size();

  auto tmp_k_arc = temp_directory_path() / "k_alt.arc";
  k_db.save(tmp_k_arc);

  REQUIRE(exists(tmp_k_arc));

  auto load_k_db = DataBase1KG{};
  load_k_db.load(tmp_k_arc);
  auto& load_k = *load_k_db.db_map.begin();

  CHECK(load_k.db_vec_snp.size() == snp_cnt);
  CHECK(load_k.db_vec_del.size() == del_cnt);
  CHECK(load_k.db_vec_ins.size() == ins_cnt);
}

TEST_CASE("Test correctness of new K"){
  using namespace Sherloc::DB;

  auto k = make_test_k_alt();

  // SNP in K
  auto snp_allele = Sherloc::SherlocMember("1", 16103, "T", "G");
  CHECK(k.find(snp_allele) == Approx(0.02f));

  // SNP not in K
  auto snp_none_allele = Sherloc::SherlocMember("1", 16103, "T", "A");
  CHECK(k.find(snp_none_allele) == Approx(0.0f));

  // INSERTION in K
  auto ins_allele = Sherloc::SherlocMember("1", 91552, "-", "T");
  CHECK(k.find(ins_allele) == Approx(0.07f));

  // INSERTION not in K
  auto ins_none_allele = Sherloc::SherlocMember("1", 9999999, "-", "AATTT");
  CHECK(k.find(ins_none_allele) == Approx(0.0f));

  // DELETION in K
  auto del_allele = Sherloc::SherlocMember("1", 83912, "AGAG", "-");
  CHECK(k.find(del_allele) == Approx(0.16f));

  // DELETION not in K
  auto del_none_allele = Sherloc::SherlocMember("1", 10000000, "AAAA", "-");
  CHECK(k.find(del_none_allele) == Approx(0.0f));
}

TEST_CASE("Test multi chr vcf file parsing"){
  auto k = make_multi_chr_k();

  auto chr1_snp_allele = Sherloc::SherlocMember("1", 55545, "C", "T");
  CHECK(k.find(chr1_snp_allele) == Approx(0.26f));

  auto chr2_snp_allele = Sherloc::SherlocMember("2", 11336, "C", "G");
  CHECK(k.find(chr2_snp_allele) == Approx(0.15f));
}