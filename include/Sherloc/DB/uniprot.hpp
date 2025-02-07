#pragma once

#include <vector>
#include <regex>
#include <string>
#include <boost/algorithm/string.hpp>
#include <map>
#include <fstream>
#include <ctime>
#include <Sherloc/DB/db.hpp>

namespace Sherloc::DB {
    
class DataBaseUniprot : public BaseDB {
    auto get_time_t(const std::string& s){
        tm t{};
        if(strptime(s.c_str(), "%Y-%m-%d", &t)){
            return mktime(&t);
        }else{
            return time_t{-1};
        }
    }
  public:
    std::map<std::string, std::string> uniprot;

    HOLMES_SERIALIZE(ar, _ver){
        ar & db_version;
        ar & db_build_time;
        ar & uniprot;
    }

    void from(const Path& filename) override {
        using namespace std::literals;
        auto mod_date_pattern = std::regex(R"token(modified="(\d+-\d{2}-\d{2})")token");
        auto var_id_pattern   = std::regex(R"token(id="(\w+)")token");

        static constexpr auto mod_date_prefix = "<entry"sv;
        static constexpr auto var_id_prefix   = "<feature"sv; // with preceding space, so can't use starswith
        static constexpr auto variant_prefix  = "<variation"sv; // with preceding space, so can't use starswith

        static constexpr auto seqvar_type     = "sequence variant"sv;

        auto is = std::ifstream{filename};
        auto latest = std::string{};
        auto latest_time = std::time_t{0};
        auto line = std::string{};
        auto matches = std::smatch{};
        auto line_num = int{0};
        while(std::getline(is, line)){
            // line that has last modified date
            if(line.starts_with(mod_date_prefix)){
                auto success = std::regex_search(line, matches, mod_date_pattern);
                if(!success){
                    SPDLOG_ERROR("<DataBaseUniprot> failed to get version on line#{}, content: ''",
                        line_num, line);
                    exit(1);
                }
                auto date = matches[1].str();
                auto mod_time = get_time_t(date);
                if(mod_time > latest_time){ // newer than latest mod time
                    latest_time = mod_time;
                    latest = std::move(date);
                }
            }

            // line that is sequence variant feature
            if(
                line.find(var_id_prefix) != std::string::npos and
                line.find(seqvar_type)   != std::string::npos
            ){
                auto success = std::regex_search(line, matches, var_id_pattern);
                if(!success){
                    SPDLOG_DEBUG("<DataBaseUniprot> variant has no id on line#{}",
                        line_num, line);
                }else{
                    auto id = matches[1].str();

                    // try to find the variation residue
                    while(std::getline(is, line)){
                        ++line_num;
                        if(line.find(variant_prefix) != std::string::npos){
                            auto residue = line.substr(line.find('>') + 1);
                            uniprot.emplace(std::move(id), residue.substr(0, residue.find('<')));
                            break;
                        }
                    }
                }
            }
            ++line_num;

            if(line_num % 2000000 == 0){
                SPDLOG_INFO("<DataBaseUniprot> processed {} lines, current db size: {}, latest mod time: {}",
                    line_num, uniprot.size(), latest);
            }
        }
        db_version = latest;

        SPDLOG_INFO("<DataBaseUniprot> Done. Total size: {}, show some entries:",
            uniprot.size());
        for(auto count = 0; auto& [id, res] : uniprot){
            SPDLOG_INFO("<DataBaseUniprot> ID: {}, residue: {}", id, res);
            ++count;
            if(count >= 10){
                break;
            }
        }
    }

    void save(const Path& filename) override {
        this->set_build_time();
        this->log_metadata("DataBaseUniprot");
        save_archive_to(*this, filename);
    }

    void load(const Path& filename) override {
        load_archive_from(*this, filename);
        this->log_metadata("DataBaseUniprot");
    }

    std::string get_codon(const std::string& str){
        auto it = uniprot.find(str);
        if(it == uniprot.end())
            return "";
        return it->second;
    }
};

}
