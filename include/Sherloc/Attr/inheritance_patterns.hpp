#pragma once

#include <map>
#include <string>
#include <Sherloc/Attr/utils.hpp>

namespace Sherloc::Attr {

struct InheritancePatterns{
    constexpr static auto patterns = make_sv_array(
        "AUTOSOMAL DOMINANT", "AUTOSOMAL RECESSIVE", "X-LINKED", "Y-LINKED"
    );
    inline static auto pattern2char(std::string_view pattern){
        static auto pattern2char_map = std::map<std::string, char>{
            {"AUTOSOMAL DOMINANT", 'D'},
            {"AUTOSOMAL RECESSIVE", 'R'},
            {"X-LINKED", 'X'},
            {"Y-LINKED", 'Y'},
            {"UNKNOWN", 'U'}
        };
        return pattern2char_map.at(std::string{pattern});
    }

    inline static auto char2abbreviation(char pattern) -> std::string {
        switch (pattern) {
            case 'D': return "AD";
            case 'R': return "AR";
            case 'X': return "XL";
            case 'Y': return "YL";
            case 'U':
            default : return "Unk";
        }
        return "Unk";
    }
};

}
