#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch/catch.hpp>

#include <fstream>
#include <vector>
#include <algorithm>
#include <string>
#include <ranges>
#include <filesystem>
#include <random>

#include <Sherloc/DB/db.hpp>
#include <Sherloc/DB/vep.hpp>
#include <Sherloc/sherloc_member.hpp>

using namespace std::filesystem;
using namespace Sherloc::DB;

TEST_CASE("VEP test"){
    // TODO: 
}