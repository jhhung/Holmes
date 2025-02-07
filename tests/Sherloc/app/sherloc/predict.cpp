#define CATCH_CONFIG_ENABLE_BENCHMARKING
#include <catch/catch.hpp>

#include <vector>
#include <algorithm>

#include <Sherloc/sherloc_member.hpp>
#include <Sherloc/Patient/patient.hpp>
#include <Sherloc/app/sherloc/sherloc_parameter.hpp>
#include <Sherloc/app/sherloc/predict.hpp>

TEST_CASE("Predict tree"){
    using namespace Sherloc::app::sherloc;
    using Sherloc::Patient::Patient;
    using Sherloc::SherlocMember;
    
    SherlocParameter para{};
    Predict predict_tree;
    Patient p;

    p.sher_mems.emplace_back("1", 1, "A", "C");
    decltype(auto) var1 = p.sher_mems.back().variants.emplace_back();

    SECTION("Go protein"){
        var1.type.emplace_back("missense_variant");
        var1.sift = true;
        var1.poly = true;
        var1.revel_score = 0.9;

        predict_tree.run(p);

        REQUIRE(var1.group.contains(122));
    }
}