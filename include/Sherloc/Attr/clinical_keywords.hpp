#pragma once

#include <map>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <Sherloc/Attr/utils.hpp>

namespace Sherloc::Attr {

struct ClinicalKeywords{
    constexpr static auto early_onset_keyword = make_sv_array(
        "EARLY", "CHILDHOOD", "NEONATAL", "INFANTILE", "FETAL");
    constexpr static auto severe_keyword = make_sv_array(
        "SEVERE", "LETHAL");
    constexpr static auto pathogenic_keyword = make_sv_array(
        "PATHOGENIC", "LIKELY_PATHOGENIC", "PATHOGENIC/LIKELY_PATHOGENIC");
    constexpr static auto benign_keyword = make_sv_array(
        "BENIGN", "LIKELY_BENIGN", "BENIGN/LIKELY_BENIGN", "BENIGN*");

    static auto contains_early_onset_keyword(std::string_view description){
        std::vector<std::string> tokens;
        boost::split(tokens, description, boost::is_any_of("_,|-"), boost::token_compress_on);
        auto onset_it = std::ranges::find(tokens, "ONSET");
        if(onset_it != tokens.end()){
            auto onset_kw_pos = std::ranges::find(early_onset_keyword, (*--onset_it).c_str());
            if(onset_kw_pos != std::end(early_onset_keyword)){
                return true;
            }
        }
        return false;
    }

    static auto contains_severe_keyword(std::string_view description){
        for(auto kw : severe_keyword){
            if(contains_str(description, kw)){
                return true;
            }
        }
        return false;
    }

    template<bool Exact = false>
    static auto contains_pathogenic_keyword(std::string_view description){
        for(auto kw : pathogenic_keyword){
            if constexpr (Exact){
                if(description == kw)
                    return true;
            } else {
                if(contains_str(description, kw))
                    return true;
            }
        }
        return false;
    }

    template<bool Exact = false>
    static auto contains_benign_keyword(std::string_view description){
        for(auto kw : benign_keyword){
            if constexpr (Exact){
                if(description == kw)
                    return true;
            } else {
                if(contains_str(description, kw))
                    return true;
            }
        }
        return false;
    }
};

}
