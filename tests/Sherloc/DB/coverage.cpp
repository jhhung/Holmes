#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch/catch.hpp>

#include <vector>
#include <algorithm>
#include <string>
#include <ranges>

#include <Sherloc/DB/coverage.hpp>

// Because this benckmark takes some time, only turn this option on when needed
// #define TURN_ON_GNOM_COVERAGE_BENCHMARK

TEST_CASE("new Gnom_coverage correctness"){
    using namespace Sherloc::app::sherloc;
    using namespace Sherloc::DB;
    using namespace std::string_literals;
    auto cov = DataBaseCoverage{};
    auto cov_vec = std::vector<Coverage>{
        Coverage{0, '0'},
        Coverage{1, '4'},
        Coverage{3, '2'},
        Coverage{70, '3'},
        Coverage{200, '4'},
        Coverage{1000, '2'},
        Coverage{2000, '3'}
    };
    cov.db_map.emplace_back(std::begin(cov_vec), std::end(cov_vec));

    auto chr = "1"s;
    
    SECTION("Basic case"){
        auto pos = size_t{5};

        CHECK(cov.find(chr, pos).status == '2');
    }

    SECTION("Edge case"){
        auto pos = size_t{70};

        CHECK(cov.find(chr, pos).status == '3');
    }

    SECTION("Head case"){
        auto pos = size_t{0};

        CHECK(cov.find(chr, pos).status == '0');
    }

    SECTION("Tail case"){
        auto pos = size_t{2001};

        CHECK(cov.find(chr, pos).status == '3');
    }
}
