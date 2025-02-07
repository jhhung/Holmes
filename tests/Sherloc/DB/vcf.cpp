#include <catch/catch.hpp>

#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <ranges>
#include <filesystem>
#include <random>

#include <Sherloc/DB/vcf.hpp>
#include <Sherloc/sherloc_member.hpp>

using namespace std::filesystem;
using namespace Sherloc::DB;

const auto child_vcf = path(DATA_PATH) / "trio/child.vcf";
const auto dad_vcf = path(DATA_PATH) / "trio/dad.vcf";
const auto mom_vcf = path(DATA_PATH) / "trio/mom.vcf";
const auto util_vcf = path(DATA_PATH) / "vcf/util.vcf";

auto count_db_size(const DataBaseVcf& db){
    auto s = size_t{};
    for(auto& [chr, map] : db.vcf_map){
        s += map.size();
    }
    return s;
}

TEST_CASE("VCF GZIPPED"){
    // TODO: I found htslib vcf parser can work with gzipped vcf file
    // here I want to check which senario might crash the parser
}

TEST_CASE("Trio VCF"){
    auto child = DataBaseVcf{};
    child.from(child_vcf);
    auto dad = DataBaseVcf{};
    dad.from(dad_vcf);
    auto mom = DataBaseVcf{};
    mom.from(mom_vcf);

    REQUIRE(count_db_size(child) == 4);
    REQUIRE(count_db_size(dad) == 4);
    REQUIRE(count_db_size(mom) == 2);

    // TODO: 
}

TEST_CASE("VCF utilities test"){
    auto vcf = HTS_VCF(util_vcf);

    auto has_value_entry = vcf.get_generic_header_value("source");
    REQUIRE(has_value_entry.has_value());
    CHECK(has_value_entry.value() == "IGSRpipeline");

    auto no_entry = vcf.get_generic_header_value("whatever");
    CHECK_FALSE(no_entry.has_value());

    auto empty_entry = vcf.get_generic_header_value("for_test");
    REQUIRE(empty_entry.has_value());
    CHECK(empty_entry.value().empty());
}

TEST_CASE("HTS VCF get_fmt"){
    SECTION("Has FORMAT"){
        auto vcf = HTS_VCF(dad_vcf);
        REQUIRE(vcf.parse_line() == HTS_VCF::OK);
        CHECK(vcf.get_fmt() == "0/1:13,13:26:99:318,0,318");
    }

    SECTION("Doesn't has FORMAT"){
        auto vcf = HTS_VCF(util_vcf);
        REQUIRE(vcf.parse_line() == HTS_VCF::OK);
        CHECK(vcf.get_fmt() == "::::");
    }
}