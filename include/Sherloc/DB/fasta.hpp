#pragma once

#include <boost/function/function_base.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <Sherloc/Attr/allele.hpp>
#include <Sherloc/DB/db.hpp>
#include <spdlog/spdlog.h>
#include <htslib/faidx.h>

namespace Sherloc::DB {

class FaidxWrapper {
private:
    faidx_t *fasta_index = nullptr;
    char* seq = nullptr;
    hts_pos_t seqlen = 0;

    inline void free_seq(){
        if(seq){
            free(seq);
        }
    }
public:
    void load_fa(const Path& fa_file){
        fai_destroy(fasta_index);
        fasta_index = fai_load(fa_file.c_str());
        if(!fasta_index){
            SPDLOG_ERROR("FaidxWrapper: Can't load faidx file, path: {}", fa_file.c_str());
            throw std::runtime_error("FaidxWrapper: Can't load faidx file");
        }
    }

    FaidxWrapper() = default;
    FaidxWrapper(const Path& fa_file){
        load_fa(fa_file);
    }

    template<Attr::IsChrType ChrType>
    auto query(const ChrType& chr, const size_t start_pos, const size_t end_pos){
        using namespace std::string_view_literals;
        free_seq();
        std::string norm_chr;
        if constexpr (std::is_same_v<ChrType, std::string>){
            norm_chr = Attr::ChrMap::norm_chr(chr);
        }else{
            norm_chr = Attr::ChrMap::idx2chr(chr);
        }
        auto query_region_str = fmt::format("{}:{}-{}",
            norm_chr, start_pos, end_pos);
        seq = fai_fetch64(fasta_index, query_region_str.c_str(), &seqlen);
        return seq ? std::string_view(seq, seqlen) : "X"sv;
    }

    ~FaidxWrapper(){
        free_seq();
        fai_destroy(fasta_index);
    }
};

class DataBaseFasta : public BaseDB {
public:
    [[gnu::deprecated("DataBaseFasta are now loaded from genome assembly, no need to build archive anymore")]]
    std::vector<std::string> fasta = std::vector<std::string>(Attr::ChrMap::approved_chr.size());

    FaidxWrapper fai_wrapper;

    HOLMES_SERIALIZE(ar, version){
        // ar & db_version;
        // ar & db_build_time;
        // ar & fasta;
    }

    [[gnu::deprecated("DataBaseFasta are now loaded from genome assembly, no need to build archive anymore")]]
    void from(const Path& filename) override {
        static constexpr auto version_prefix = std::string_view{"chromosome:"};
        auto is = std::ifstream{filename};
        auto seq = std::string{};
        size_t current_chr_idx = 0;
        bool skip_chr = false;

        while(std::getline(is, seq )) {
            if(seq.starts_with('#'))
                continue;
            if(seq.starts_with('>')) {
                auto chr_str = seq.substr(1);
                auto chr = chr_str.substr(0, chr_str.find(' '));
                
                // set version
                if(db_version.empty()){
                    auto ver_pos = chr_str.find(version_prefix);
                    if(ver_pos != std::string::npos){
                        chr_str = chr_str.substr(ver_pos + version_prefix.size());
                        db_version = chr_str.substr(0, chr_str.find(':'));
                    }
                }
                try{
                    current_chr_idx = Attr::ChrMap::chr2idx(chr);
                    skip_chr = false;
                }catch(std::out_of_range& e){
                    SPDLOG_WARN("Skip chromosome '{}'", chr);
                    skip_chr = true;
                }
                continue;
            }
            if(skip_chr){
                continue;
            }
            fasta[current_chr_idx] += seq;
        }

        for(auto idx = 0; auto& seq : fasta){
            std::cout << fmt::format("Seqname: {}, size: {}\n",
                Attr::ChrMap::idx2chr(idx), seq.size());
            ++idx;
        }
    }

    [[gnu::deprecated("DataBaseFasta are now loaded from genome assembly, no need to build archive anymore")]]
    void save(const Path& filename) override {
        this->set_build_time();
        this->log_metadata("DataBaseFasta");
        save_archive_to(*this, filename);
    }

    void load(const Path& filename) override {
        // load faidx file
        fai_wrapper.load_fa(filename);

        // // [DEPRECATED] load boost binary archive
        // load_archive_from(*this, filename);
        // this->log_metadata("DataBaseFasta");
    }

    template<Attr::IsChrType ChrType>
    [[deprecated("now the reference seq can be accessed using fai_wrapper::query")]]
    char get_ref(const ChrType& chr, size_t pos) {
        if constexpr (std::is_same_v<ChrType, std::string>) {
            return get_ref(Attr::ChrMap::chr2idx(chr), pos);
        } else if constexpr (std::is_same_v<ChrType, size_t>) {
            return fasta.at(chr).at(pos-1);
        }
        return '\0'; // TODO: never get here anyway, (in c++23, we can use std::unreachable())
    }

    inline bool check(const Attr::IsChrType auto& chr, size_t start, size_t end, const std::string& ref) {
        return fai_wrapper.query(chr, start, end) == ref;
    }

    inline bool check_base(const Attr::IsChrType auto& chr, size_t pos, char ref) {
        return fai_wrapper.query(chr, pos, pos)[0] == ref;
    }
};

}
