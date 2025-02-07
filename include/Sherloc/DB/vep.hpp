#pragma once

#include <ranges>
#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <optional>
#include <boost/algorithm/string/split.hpp>
#include <Sherloc/Attr/utils.hpp>
#include <Sherloc/Attr/allele.hpp>
#include <Sherloc/variant.hpp>
#include <Sherloc/DB/vcf.hpp>
#include <Sherloc/DB/db.hpp>
#include <Sherloc/sherloc_member.hpp>
#include <spdlog/spdlog.h>
#include <nlohmann/json.hpp>

namespace Sherloc::DB {

using Path = std::filesystem::path;

class VEP : BaseDB {
public:
    struct Table;

    using PosType = u_int32_t;
    using IndexType = size_t;
    using CacheType = std::map<size_t, Table>;

    struct Table{
        std::vector<PosType> positions;
        std::vector<std::string> refs;
        std::vector<std::string> alts;
        std::vector<std::string> records;

        HOLMES_SERIALIZE(ar, version){
            ar & positions;
            ar & refs;
            ar & alts;
            ar & records;
        }

        Table() = default;

        Table(int reserved_size) {
            positions.reserve(reserved_size);
            refs.reserve(reserved_size);
            alts.reserve(reserved_size);
            records.reserve(reserved_size);
        }

        template<class Str>
        void add_row(const HTS_VCF::VCF_Record& allele, Str&& record){
            positions.emplace_back(allele.pos);
            refs.emplace_back(allele.ref);
            alts.emplace_back(allele.alt);

            records.emplace_back(std::forward<Str>(record));
        }

        [[nodiscard]] auto is_sorted() const {
            return std::ranges::is_sorted(positions);
        }
    };

    struct CacheMeta{
        std::map<std::string, std::string> meta;
        HOLMES_SERIALIZE(ar, version){
            ar & meta;
        }
    };
private:
    Attr::HeaderIndexType header_index;
    CacheMeta cache_meta;
    Path out_dir;
    Table current_table;
    std::string current_chr = "";
    static constexpr int reserve_size = 100000;
    static constexpr std::string_view meta_filename = "meta.arc";

    auto try_load_chr(const std::string& normed_chr){
        // assume the chr is normalized
        auto target_chr_file = out_dir / fmt::format("{}.arc", normed_chr);
        if(!std::filesystem::exists(target_chr_file)){
            return false;
        }
        load_archive_from(current_table, target_chr_file);
        return true;
    }
public:
    HOLMES_SERIALIZE(ar, version){
        ar & db_version;
        ar & db_build_time;
        ar & header_index;
        ar & cache_meta;
    }

    /**
    * @brief Reads the header of a VEP tab-delimited output file and returns a map of column names to their indices.
    * 
    * @param is The input stream to read from.
    * @return std::map<std::string, size_t> A map of column names to their indices.
    * 
    * This function is used for parsing the header of a VEP tab-delimited output file. It reads the header line and
    * extracts the column names, then returns a map of column names to their indices.
    */
    static auto read_header(std::istream& is){
        static constexpr auto tab_delimiter = Attr::delimiter('\t');

        std::string buffer;

        // parse header
        Attr::HeaderIndexType vep_header_index;
        while(is.peek() == '#'){
            is.get();
            if(is.peek() == '#') // this line (something like ## HGVSp : ...) is not the header line
                std::getline( is, buffer );
            else{ // this line (something like #Uploaded_variation ....) is the header line
                std::getline( is, buffer );
                vep_header_index = Attr::make_header_index(buffer, tab_delimiter);
            }
        }
        return vep_header_index;
    }

    VEP() = default;

    void set_output_dir(const Path& dir){
        out_dir = dir;
    }

    template <class Container>
    static auto parse_vcf_into(
        HTS_VCF& vcf,
        Container& container,
        int reserved_size = 0 /* this variable is not used if not parsing into cache */
    ){
        constexpr auto pipe_delimiter = Attr::delimiter('|');
        HTS_VCF::VCF_Status vcf_status;

        // parse vep vcf header
        auto desc = vcf.get_header_info_description("CSQ");
        SPDLOG_INFO("VEP VCF CSQ description: {}", desc);
        auto start_pos = desc.find("Format: ") + 8;
        auto end_pos = desc.find('"', start_pos);
        Attr::HeaderIndexType vep_header_index = Attr::make_header_index(
            desc.substr(start_pos, end_pos - start_pos), pipe_delimiter);

        // add records
        int line = 0;
        while((vcf_status = vcf.parse_line()) != HTS_VCF::VCF_Status::VCF_EOF){
            switch (vcf_status) {
                case HTS_VCF::VCF_Status::OK:
                    break;
                case HTS_VCF::VCF_Status::RECORD_NO_ALT:
                    // Assume not to happen
                    continue;
                default:
                    SPDLOG_ERROR("HTS_VCF parsing error, status code: {}", int(vcf_status));
                    exit(1);
            }

            if constexpr (std::is_same_v<Container, CacheType>){ // parsing vcf into vep cache
                auto [it, success] = container.emplace(
                    Attr::ChrMap::chr2idx(vcf.record.chr),
                    reserved_size);

                it->second.add_row(
                    vcf.record,
                    vcf.info_str("CSQ")
                       .value_or("")
                );
            } else { // parsing into sherloc_member vector
                container.at(std::stoul(vcf.get_ID())).variants = 
                    Variant::make_variants(
                        vep_header_index,
                        vcf.info_str("CSQ")
                           .value_or("") // views::split ranges will have 0 size given an empty string
                    );
            }

            ++line;
            if(line % 100000 == 0){
                SPDLOG_INFO("VEP parsed {} lines.", line);
            }
        }
        return vep_header_index;
    }

    void from(const Path& filename) override {
        HTS_VCF vep_vcf{filename, true, false, true, false};
        this->db_version = vep_vcf
            .get_generic_header_value("VEP")
            .value_or("None");
        CacheType cache;

        header_index = parse_vcf_into(vep_vcf, cache, reserve_size);

        // check all chromosome are sorted by position
        for(auto& [chr, table] : cache){
            auto chr_str = Attr::ChrMap::idx2chr(chr);
            if(!table.is_sorted()){
                SPDLOG_ERROR(
                    "chr{} is not sorted! pleaset sort the vcf first.",
                    chr_str);
                exit(1);
            }else{
                SPDLOG_DEBUG("chr{} is sorted.", chr_str);
            }
        }

        // save to out dir
        if(!std::filesystem::exists(out_dir)){
            std::filesystem::create_directories(out_dir);
        }
        for(auto& [chr, table] : cache){
            auto chr_str = Attr::ChrMap::idx2chr(chr);
            save_archive_to(table, out_dir / fmt::format("{}.arc", chr_str));
        }
    }

    void load(const Path& dirname) override {
        set_output_dir(dirname);
        load_archive_from(*this, dirname / meta_filename);
        this->log_metadata("VEPCache");
        current_chr = ""; // clean the current chr
    }

    void save(const Path& dirname) override {
        if(!std::filesystem::exists(dirname)){
            std::filesystem::create_directories(dirname);
        }
        this->set_build_time();
        this->log_metadata("VEPCache");
        save_archive_to(*this, dirname / meta_filename);
    }

    [[nodiscard]] auto find(const SherlocMember& sher_mem) 
        -> std::optional<std::vector<std::string>::const_iterator>
    {
        if(sher_mem.chr != current_chr){
            if(!try_load_chr(sher_mem.chr)){ // can't not load the cache chr
                return std::nullopt;
            }
            current_chr = sher_mem.chr;
        }

        auto [s_it, e_it] = std::ranges::equal_range(current_table.positions, sher_mem.pos);
        if(s_it == e_it){
            return std::nullopt;
        }

        auto s_idx = std::distance(current_table.positions.begin(), s_it);
        auto e_idx = std::distance(current_table.positions.begin(), e_it);

        for(auto idx = s_idx; idx != e_idx; ++idx){
            if(current_table.refs[idx] == sher_mem.ref and current_table.alts[idx] == sher_mem.alt){
                return current_table.records.cbegin() + idx;
            }
        }
        return std::nullopt;
    }
    
    /**
     * @brief Tries to find if a SherlocMember allele is in the cache. 
     * If found, they will be inserted into SherlocMember and return true, else return false.
     * 
     * @param sher_mem The SherlocMember object to insert the variant into.
     * @return true If the allele is found in the cache and inserted into SherlocMember.
     * @return false If the allele is not found in the cache.
     */
    auto try_insert_into(SherlocMember& sher_mem) {
        auto it = find(sher_mem);
        if(it.has_value()){
            sher_mem.variants = Variant::make_variants(header_index, *(it.value()));
            return true;
        }
        return false;
    }
};

class VEPRunner {
private:
    std::vector<std::string> options = {};
    std::map<std::string, std::string> plugins = {};
    
    Path vep_repo_dir;
    Path vep_data_dir;
    Path vep_executable;
    Path filter_vep_executable;
    Path assembly_file;
public:
    VEPRunner() = default;

    VEPRunner(const Path& config_file) {
        auto config = Attr::load_json(config_file,
            fmt::format("vep runner: config_file '{}' can't be opened!",
                config_file.c_str()));

        vep_repo_dir = config["vep_repo_dir"].get<std::string>();
        if(vep_repo_dir.empty()){
            vep_repo_dir = Path(std::getenv("HOME")) / "ensembl-vep";
            SPDLOG_WARN("vep_repo_dir: is empty. Assume to be `$HOME/ensembl-vep` = {}",
                vep_repo_dir.c_str());
            Attr::check_exist(vep_repo_dir, fmt::format("'{}'", vep_repo_dir.c_str()));
        }

        vep_data_dir = config["vep_data_dir"].get<std::string>();
        if(vep_data_dir.empty()){
            vep_data_dir = Path(std::getenv("HOME")) / ".vep/for_vep";
            SPDLOG_WARN("vep_data_dir: is empty. Assume to be `$HOME/.vep/for_vep` = {}",
                vep_data_dir.c_str());
            Attr::check_exist(vep_data_dir, fmt::format("'{}'", vep_data_dir.c_str()));
        }

        vep_executable = vep_repo_dir / "vep";
        filter_vep_executable = vep_repo_dir / "filter_vep";

        assembly_file = config["assembly"].get<std::string>();
        if(assembly_file.is_relative())
            assembly_file = vep_repo_dir / assembly_file;

        Attr::check_exist(vep_executable, fmt::format("'{}'", vep_executable.c_str()));
        Attr::check_exist(filter_vep_executable, fmt::format("'{}'", filter_vep_executable.c_str()));
        Attr::check_exist(assembly_file, fmt::format("'{}'", assembly_file.c_str()));

        options = config["options"].get<std::vector<std::string>>();
        
        // make plugins cmd
        for(auto&& [plugin_name, plugin_files] : 
            config["plugins"].get<std::map<std::string, std::vector<std::string>>>())
        {
            auto plugin_cmd = fmt::format("--plugin {}", plugin_name);
            if(!plugin_files.empty()){ // plugin needs input file(s)
                plugin_cmd += fmt::format(
                    ",{}", fmt::join(plugin_files | std::views::transform([&vep_data_dir = vep_data_dir, &plugin_name=plugin_name](const auto& file){
                        // FIXME: alphamissense special case
                        if (plugin_name == "AlphaMissense"){
                            return file;
                        }
                        auto file_p = Path(file);
                        return (file_p.is_relative() ? vep_data_dir / file_p : file_p).string();
                    }), ","));
            }
            plugins[plugin_name] = plugin_cmd;
        }
    }

    [[nodiscard]] auto make_cmd(
        const Path& vep_inputfile,
        const Path& vep_outputfile,
        const std::vector<std::string>& gene_list,
        int thread_num,
        bool is_grch37 = false
    ) const {
        auto vep_args = options;
        auto run_filter = gene_list.size() > 0;

        // input / output / assembly file / thread
        // assembly file is either absolute path or relative path under vep_repo_dir
        vep_args.emplace_back(
            fmt::format("-i {} -o {} --assembly {} --fasta {} --fork {}",
                vep_inputfile.c_str(),
                run_filter ? "STDOUT" : vep_outputfile.c_str(),
                is_grch37 ? "GRCh37" : "GRCh38",
                assembly_file.c_str(),
                std::max(thread_num, 1) // positive value
            )
        );

        // push vep plugins command
        for(auto&& [plugin, cmd] : plugins) {
            SPDLOG_INFO("[run vep] VEP uses {} plugin", plugin);
            vep_args.emplace_back(cmd);
        }

        SPDLOG_INFO("[run vep] filter_vep: {}", (run_filter ? "on" : "off"));
        auto filter_vep_cmd = std::string{};
        if(run_filter){
            filter_vep_cmd = fmt::format(
                "| {} -o {} --filter \"SYMBOL {} {}\" --force_overwrite",
                    filter_vep_executable.c_str(),
                    vep_outputfile.c_str(),
                    (gene_list.size() == 1 ? "is" : "in"),
                    fmt::join(gene_list, ",")
            );
        }
        auto vep_cmd = fmt::format(
            "{} {} {}",
                vep_executable.c_str(),
                fmt::join(vep_args, " "),
                filter_vep_cmd
        );
        return vep_cmd;
    }

    [[nodiscard]] auto get_assembly_file() const {
        return assembly_file;
    }
};

}
