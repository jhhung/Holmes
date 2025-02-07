#include <catch/catch.hpp>

#include "Sherloc/Attr/allele.hpp"

TEST_CASE("Test Allele"){
    // TODO: 
}

TEST_CASE("Test ChrMap chr2idx"){
    using namespace Sherloc::Attr;

    CHECK(ChrMap::chr2idx("chr1") == 0);
    CHECK(ChrMap::chr2idx("chrX") == 22);
    CHECK(ChrMap::chr2idx("chrY") == 23);

    CHECK(ChrMap::chr2idx("1") == 0);
    CHECK(ChrMap::chr2idx("3") == 2);
    CHECK(ChrMap::chr2idx("X") == 22);
    CHECK(ChrMap::chr2idx("Y") == 23);

    CHECK_THROWS(ChrMap::chr2idx("chrA"));
    CHECK_THROWS(ChrMap::chr2idx("ch1"));
    CHECK_THROWS(ChrMap::chr2idx("chr23"));
    CHECK_THROWS(ChrMap::chr2idx("chr100"));
    CHECK_THROWS(ChrMap::chr2idx("100"));
    CHECK_THROWS(ChrMap::chr2idx(""));
}

TEST_CASE("Test ChrMap idx2chr"){
    using namespace Sherloc::Attr;

    CHECK(ChrMap::idx2chr(1) == "2");
    CHECK(ChrMap::idx2chr(0) == "1");
    CHECK(ChrMap::idx2chr(23) == "Y");

    CHECK_THROWS(ChrMap::idx2chr(24));
    CHECK_THROWS(ChrMap::idx2chr(100));
    CHECK_THROWS(ChrMap::idx2chr(-1));
}