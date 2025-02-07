#pragma once

#include <string>
#include <map>
#include <vector>
#include <algorithm>
#include <Sherloc/Attr/utils.hpp>

namespace Sherloc::Attr {

struct Allele {
    // TODO: I want to integrate all "allele" related data
    // like [chr, pos, ref, alt] tuple to one class, and unify 
    // vcf-style and vep-style representation
    size_t chr_idx;
    size_t pos;
    std::string ref;
    std::string alt;

    enum DeNovo: int {
        Unknown,
        IsDeNovo,
        NotDeNovo
    };


    enum ViewType : int {
        VCF,
        VEP
    };
    
    template<ViewType VT = ViewType::VCF>
    struct View {
        size_t chr_idx;
        size_t pos;
        std::string_view ref;
        std::string_view alt;

        View() = default;
        View(const View &) = default;
        View(View &&) = default;
        View &operator=(const View &) = default;
        View &operator=(View &&) = default;

        View(const Allele &allele)
            : chr_idx(allele.chr_idx), pos(allele.pos), ref(allele.ref), alt(allele.alt) {
            if constexpr (VT == ViewType::VEP){
                static constexpr auto null_base = std::string_view{"-"};
                if (ref.size() != alt.size()) {
                    auto same_prefix_len = std::distance(
                        std::begin(ref), std::ranges::mismatch(ref, alt).in1);
                    ref.remove_prefix(same_prefix_len);
                    alt.remove_prefix(same_prefix_len);
                    if (ref.empty())
                    ref = null_base;
                    if (alt.empty())
                    alt = null_base;
                }
            }
        }

        auto operator<=>(const auto& rhs) const {
            return 
                std::tie(chr_idx, pos, ref, alt) <=>
                std::tie(rhs.chr_idx, rhs.pos, rhs.ref, rhs.alt);
        }
    };

    auto operator<=>(const auto& rhs) const {
        return 
            std::tie(chr_idx, pos, ref, alt) <=>
            std::tie(rhs.chr_idx, rhs.pos, rhs.ref, rhs.alt);
    }
};

template<typename T>
concept IsChrType = IsAnyOf<T, std::string, size_t>;

struct ChrMap {
private:
    static auto _make_chr_map() {
        std::map<std::string, size_t, std::less<void>> chr_map{};
        for(size_t idx = 0; auto chr : approved_chr){
            chr_map.emplace(chr, idx);
            chr_map.emplace("chr" + std::string{chr}, idx);
            ++idx;
        }
        return chr_map;
    }
public:

    constexpr static auto approved_chr = make_sv_array(
        "1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
        "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
        "21", "22", "X", "Y");

    static auto chr2idx(const std::string& chr) {
        static const auto chr_map{_make_chr_map()};
        return chr_map.at(chr);
    }

    constexpr static auto idx2chr(size_t idx) {
        return approved_chr.at(idx);
    }

    static auto norm_chr(const std::string& chr) {
        return idx2chr(chr2idx(chr));
    }
};

}
