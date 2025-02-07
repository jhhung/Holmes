#pragma once

#include <string>
#include <array>
#include <map>
#include <filesystem>
#include <fstream>
#include <optional>
#include <spdlog/spdlog.h>
#include <boost/algorithm/string/split.hpp>
#include <nlohmann/json.hpp>

namespace Sherloc::Attr {

using Path = std::filesystem::path;
using HeaderIndexType = std::map<std::string, size_t, std::less<>>;

template <typename T, typename... U>
concept IsAnyOf = (std::same_as<T, U> || ...);

/**
 * @brief Creates a std::array of std::string_view from a variadic list of arguments.
 * the size of array can be auto deduced
 * 
 * @tparam SVs Variadic template parameter pack of std::string_view types.
 * @param svs Variadic list of std::string_view arguments.
 * @return std::array<std::string_view, sizeof...(svs)> An std::array of std::string_view.
 */
template<typename ... SVs>
constexpr decltype(auto) make_sv_array(SVs&& ... svs){
    return std::array<std::string_view, sizeof...(svs)>{std::forward<SVs>(svs)...};
}

/**
 * @brief Check if a string contains another string.
 * 
 * This function checks if the reference string contains the query string.
 * 
 * @param ref_str The reference string to search in.
 * @param query_str The query string to search for.
 * @return true if the reference string contains the query string, false otherwise.
 * 
 * @note This function is case-sensitive.
 */
inline auto contains_str(std::string_view ref_str, std::string_view query_str){
    return ref_str.find(query_str) != std::string::npos;
}

/**
 * @brief Convert a character to uppercase. Can be used as std::ranges::for_each parameter to 
 * convert string to upper case inplace.
 * 
 * This function converts the given character to uppercase using the `toupper` function.
 * 
 * @param c The character to convert to uppercase.
 */
inline void as_upper(char& c){
    c = std::toupper(c);
}

/**
 * @brief A factory function of single character comparison lambda
 * I write this function because boost::is_any_of takes a set of chars
 * but sometime we only need one char
 * 
 * @param d Delimiter
 * @return auto Lambda function that return if a char is d or not
 */
constexpr inline auto delimiter(char d){
    return [d](auto c) constexpr -> bool { return c == d; };
}

template<typename T>
concept Findable = 
    requires(T x, std::string str) {
        { str.find(x) } -> std::convertible_to<decltype(std::string::npos)>;
    };

/**
 * @brief Split a string into two parts based on a delimiter.
 *
 * This function takes a string `str` and a delimiter `delimiter` and returns
 * a pair of strings representing the two parts of the string split by the delimiter.
 * The type of the delimiter can be either a string-like type (such as `const char*`
 * or `std::string`) or a single character.
 *
 * @tparam Delimiter The type of the delimiter. This can be either a string-like type
 * (such as `const char*` or `std::string`) or a single character.
 *
 * @param str The input string to split.
 * @param delimiter The delimiter to split the string by.
 *
 * @return A pair of strings representing the two parts of the input string split by
 * `delimiter`.
 */
template <Findable Delimiter>
inline auto explode(const std::string& str, const Delimiter& delimiter){
    auto split_pos = str.find(delimiter);
    size_t len;
    if constexpr (std::is_convertible_v<Delimiter, std::string_view>) { // string_like type (const char*, string, etc..)
        len = std::string_view(delimiter).size();
    } else { // char
        len = 1;
    }
    return split_pos != std::string::npos ?
        std::make_pair(str.substr(0, split_pos), str.substr(split_pos + len)) :
        std::make_pair(str                     , "");
};

/**
 * @brief Parse VEP cds/amino acid position format.
 *
 * This function takes a string representing the position of a variant in a gene
 * in VEP cds/amino acid position format and returns an unsigned long integer
 * representing the position. The input string can have one of two formats:
 * a single number or two numbers separated by a hyphen. If the input string has
 * the format of a single number, the function returns the unsigned long integer
 * value of the number. If the input string has the format of two numbers separated
 * by a hyphen, the function returns the unsigned long integer value of the second
 * number if the first number is "?", or the unsigned long integer value of the first
 * number otherwise.
 *
 * @param position A string representing the position of a variant in a gene in
 * VEP cds/amino acid position format.
 *
 * @return An unsigned long integer representing the position of the variant in
 * the gene.
 */
inline auto parse_vep_pos(const std::string& position){
    // vep cds/aa pos might be: 
    // 1) one number
    // 2) 2 numbers sep by '-'
    //  2A) first is ?
    //  2B) second is ?
    auto&& [pos1, pos2] = explode(position, '-');
    if(pos1 != "?"){
        return std::stoul(pos1);
    }
    return std::stoul(pos2);
}

inline void check_exist(const Path& p, const std::string& fallback_msg){
    if(!std::filesystem::exists(p)){
        SPDLOG_ERROR("File '{}' not exist!",
            p.c_str());
        throw std::runtime_error(fallback_msg);
    }
}

template<typename T>
concept BoostDelimiter = 
    requires(T x, char c) {
        { x(c) } -> std::convertible_to<bool>;
    };

template <BoostDelimiter Delimiter>
inline auto make_header_index(std::string_view header_line, Delimiter&& del){
    auto index = HeaderIndexType{};
    auto cols = std::vector<std::string>{};
    boost::split(cols, header_line, std::forward<Delimiter>(del));
    for(auto idx = size_t{0}; auto& col : cols){
        index.emplace(col, idx);
        SPDLOG_DEBUG("col:idx = {}:{}", col, idx);
        ++idx;
    }
    return index;
}

inline auto load_json(
    const Path& json_file,
    const std::string& fallback = "",
    bool throw_when_failed = false
){
    nlohmann::json js;
    {
        auto ifs = std::ifstream{json_file};
        if(!ifs.is_open()){
            SPDLOG_ERROR(fallback);
            if(throw_when_failed){
                throw std::runtime_error(fmt::format("can't open {}", json_file.c_str()));
            }else{
                exit(1);
            }
        }
        ifs >> js;
    }
    return js;
}

}