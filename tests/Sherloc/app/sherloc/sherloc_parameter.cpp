#include <catch/catch.hpp>
#include <Sherloc/app/sherloc/sherloc_parameter.hpp>

TEST_CASE("Sherloc parameters initialization"){
    using namespace Sherloc::app::sherloc;

    decltype(auto) para = Sherloc::app::sherloc::SherlocParameter::get_paras();

    CHECK(para.coverage_high == 8);
    CHECK(para.coverage_low == 1);

    CHECK(para.score_table[93] == -5.);
    CHECK(para.score_table[135] == 1.);

    CHECK(para.score_table[17] == 4.);
    CHECK(para.score_table[21] == -5.);
    CHECK(para.score_table[44] == 1.5);

    CHECK(para.score_table[211] == 0.);
    CHECK(para.score_table[205] == 4.);

    CHECK(para.score_table[23] == 2.5);
    CHECK(para.score_table[188] == 0.5);
}