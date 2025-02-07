#include <catch/catch.hpp>

#include <string>
#include <Sherloc/DB/fasta.hpp>
#include <spdlog/fmt/fmt.h>

// TEST_CASE("Test htslib faidx query") {
//     using namespace std::literals;
//
//     // TODO: this is so simple and intuitive, consider refactor the DatabaseFasta with htslib faidx
//     // instead of storing the uncompressed fasta as boost archive
//
//     auto fa = "<path/to/some.fasta.gz>"s;
//     faidx_t *fasta_index;
//     char *seq;
//     int seq_len;

//     // Load the FASTA file
//     fasta_index = fai_load(fa.c_str());
//     REQUIRE(fasta_index != nullptr);

//     // Fetch a sequence
//     seq = fai_fetch(fasta_index, "1:1000-2000", &seq_len);
//     REQUIRE(seq != nullptr);

//     // Print the sequence
//     fmt::print("1:1000-2000 = '{}'\n", seq);

//     // Clean up
//     free(seq);
//     fai_destroy(fasta_index);
// }