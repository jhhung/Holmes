#pragma once

#include <utility>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <Sherloc/DB/db.hpp>
#include <Sherloc/DB/fasta.hpp>
#include <Sherloc/variant.hpp>
#include <spdlog/spdlog.h>
#include <optional>

namespace Sherloc::DB {

class Transcript
{
  public:

    size_t start = 0;

    size_t end = 0;

    std::string id = "";

    char strand = '+';

    std::vector< std::pair< int, int > > exon;

    template< class Archive >
    void serialize( Archive& ar, const unsigned int version )
	{
        ar & start;
        ar & end;
        ar & id;
        ar & strand;
        ar & exon;
    }

    Transcript() = default;

    Transcript( const std::vector< std::string >& vec ): start( std::stoi(vec[3]) ), end( std::stoi(vec[4]) )
    {
        id = Variant::feature_normalize(vec[11]);
        strand = vec[6][0];
        if(strand != '+' and strand != '-')
            throw std::runtime_error(fmt::format("Transcript, strange strand {}", strand));
    }

    void set_exon( const std::vector< std::string >& vec )
    {
        exon.emplace_back( std::stoi(vec[3]), std::stoi(vec[4]) );
    }
};

/**
 * @brief Use Gencode GTF "Basic gene annotation  CHR" file
 *  File name should be like: gencode.vXXX.basic.annotation.gtf
 *  Field: 
 *  idx[0] 1 	chromosome name 	chr{1,2,...,22,X,Y,M} or GRC accession
 *  idx[1] 2 	annotation source 	{ENSEMBL,HAVANA}
 *  idx[2] 3 	feature type 	{gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}
 *  idx[3] 4 	genomic start location 	integer-value (1-based)
 *  idx[4] 5 	genomic end location 	integer-value
 *  idx[5] 6 	score(not used) 	.
 *  idx[6] 7 	genomic strand 	{+,-}
 *  idx[7] 8 	genomic phase (for CDS features) 	{0,1,2,.}
 *  idx[8] 9 	additional information as key-value pairs
 *      TODO: Maybe we need a more general way to parse GTF attribute field?
 *            link: https://www.gencodegenes.org/pages/data_format.html
 *      idx[9] value of 'gene_id'
 *      idx[11] value of 'transcript_id' if feature type is transcript/exon
 */
class DataBaseGTF : public BaseDB {
  public:
    std::map<std::string, Transcript> txp_map;

    enum class GTFSource {
        EnsemBl,
        RefSeq
    };

    HOLMES_SERIALIZE(ar, _ver){
        ar & db_version;
        ar & db_build_time;
        ar & txp_map;
    }

    void parse_gtf(const Path& gtf_filename, GTFSource source){
        std::ifstream is(gtf_filename);
        std::string str;
        int id = 0;
        while(std::getline(is, str)) {
            if( str[0] == '#' ) {
                switch (source) {
                    using enum GTFSource;
                    case EnsemBl:
                        {
                            static constexpr std::string_view key{"genebuild-last-updated"};
                            auto pos = str.find(key);
                            if(pos != std::string::npos){
                                db_version += fmt::format("Ensembl: [{}], ",
                                    str.substr(pos + key.size() + 1));
                            }
                        }
                        break;
                    case RefSeq:
                        {
                            static constexpr std::string_view key{"annotation-source"};
                            auto pos = str.find(key);
                            if(pos != std::string::npos){
                                db_version += fmt::format("RefSeq: [{}], ",
                                    str.substr(pos + key.size() + 1));
                            }
                        }
                        break;
                }
                continue;
            }
            id++;
            if(id % 100000 == 0){
                SPDLOG_INFO("processed {} lines...", id);
            }
            std::vector< std::string > vec;
            boost::split( vec, str, boost::is_any_of( "\t ;\"" ), boost::token_compress_on );
            auto ens_txp_id = Variant::feature_normalize(vec[11]);
            if(vec[2] == "transcript")
                txp_map.emplace(ens_txp_id, vec);

            if(vec[2] == "exon")
                txp_map.find(ens_txp_id)->second.set_exon(vec);
        }
    }

    void from(const Path& filename) override {
        auto gtf_input_config = Attr::load_json(filename);

        db_version = "";
        SPDLOG_INFO("Parsing Ensembl GTF file...");
        parse_gtf(gtf_input_config["ensembl"].get<std::string>(), GTFSource::EnsemBl);
        SPDLOG_INFO("Parsing RefSeq GTF file...");
        parse_gtf(gtf_input_config["refseq"].get<std::string>(), GTFSource::RefSeq);
    }

    void save(const Path& filename) override {
        this->set_build_time();
        this->log_metadata("DataBaseGTF");
        save_archive_to(*this, filename);
    }

    void load(const Path& filename) override {
        load_archive_from(*this, filename);
        this->log_metadata("DataBaseGTF");
    }

    /**
     * @brief Find Transcript by Ensembl Gene and Transcipt ID
     * 
     * @param gene_id DEPRECATED, placeholder
     * @param trans_id 
     * @return std::optional<std::reference_wrapper<const Transcript>> 
     */
    inline std::optional<std::reference_wrapper<const Transcript>>
    find_transcript(const std::string& gene_id, const std::string& trans_id){
        auto it = txp_map.find(trans_id);
        if(it == txp_map.end())
            return std::nullopt;
        return std::cref(it->second);
    }

    std::pair< int, int > get_last_exon( const std::string& gene_id, const std::string& trans_id )
    {
        auto trans_opt = find_transcript(gene_id, trans_id);
        if(!trans_opt)
            return {0, 0};
        decltype(auto) trans = trans_opt.value().get();
        return trans.strand == '+' ? trans.exon.back() : trans.exon.front();
    }

    /**
     * @brief DEPRECATED: We found a VEP plugin that can do this, and with more robust criteria
     * Checks if a variant that is the target of nonsense-mediated decay 
     * (NMD) has chance to **ESCAPE** this mechanism.
     * 
     * This function checks if a variant at a given position in a transcript which is the 
     * target of nonsense-mediated decay (which is a mechanism that degrades transcripts with
     * premature termination codons (PTCs) to prevent the production of truncated proteins), 
     * have any chance to ESCAPE this mechanism.
     * 
     * @param gene_id The Ensembl gene ID of the gene containing the transcript.
     * @param trans_id The Ensembl transcript ID of the transcript.
     * @param pos The 0-based position of the variant in the transcript.
     * @return true if the variant might escape NMD, false otherwise.
     */
    [[gnu::deprecated]] bool might_escape_nmd( const std::string& gene_id, const std::string& trans_id, const size_t& pos )
    {
        auto trans_opt = find_transcript(gene_id, trans_id);
        if(!trans_opt.has_value())
            return false;
        decltype(auto) trans = trans_opt.value().get();
        auto is_in_this_exon = [pos](const std::pair<int, int>& exon){
            return pos > exon.first and pos < exon.second;
        };
        auto exon_it = std::ranges::find_if(trans.exon, is_in_this_exon);
        // if not found, than maybe this position has some problem
        // this shouldn't happen 'cause vep annotate this variant as "stop_gained" or "frameshift_variant"
        if(exon_it == trans.exon.end()){
            SPDLOG_DEBUG("GTF::check_nmd(): pos is not in any exon");
            return false;
        }
        if(trans.exon.size() <= 1){
            SPDLOG_DEBUG("GTF::check_nmd(): transcript has {} exons", trans.exon.size());
            return true; // only one exon means this is the last exon
        }
        int exon_idx = std::distance(trans.exon.begin(), exon_it);

        // the following condition aims to check if a null variant might escape NMD
        //  we use last 15 codon (45 bases) as the preultimate threshold
        if(trans.strand == '+'){
            if( ( exon_idx == trans.exon.size()-1 )
                or
                ( exon_idx == trans.exon.size()-2 and trans.exon[exon_idx].second - pos < 45)) 
                return true;
        }else{
            if( ( exon_idx == 0 )
                or
                ( exon_idx == 1 and pos - trans.exon[exon_idx].first < 45)) 
                return true;
        }
        return false;
    }

    /**
     * @brief Checks if a variant affects the highly conserved acceptor AG sequence in a transcript.
     * 
     * This function checks if a variant at a given position in a transcript affects the highly
     * conserved acceptor AG sequence, which is located at the 3' end of introns and is essential for
     * proper splicing of pre-mRNA. The function only checks variants that are within 2 bases of
     * the acceptor AG sequence.
     * 
     * @param gene_id The Ensembl gene ID of the gene containing the transcript.
     * @param trans_id The Ensembl transcript ID of the transcript.
     * @param chr The chromosome where the variant is located.
     * @param pos The 0-based position of the variant in the transcript.
     * @param fa A reference to a `Fasta` object that provides access to the reference genome.
     * @return true if the variant affects the acceptor AG sequence, false otherwise.
     */
    bool check_acceptor( 
        const std::string& gene_id, const std::string& trans_id,
        const std::string& chr, const size_t& pos, DataBaseFasta& fa, Variant& variant)
    {
        auto trans_opt = find_transcript(gene_id, trans_id);
        if(!trans_opt){
            // Maybe the hgvs is in RefSeq format
            // Use quick workaround by directly finding '-2A' or '-1G' in HGVSc
            return 
                (variant.hgvsc.find("-2A") != std::string::npos) or
                (variant.hgvsc.find("-1G") != std::string::npos);
        }
        decltype(auto) trans = trans_opt.value().get();

        // Check if the variant position is in the highly conserved acceptor AG 
        // |----intron-----AG|******exon****| ... 
        //                 ^^ (these two position)
        // TODO: But!!!! 
        // rule EV0198:
        //      Variant in donor GT or acceptor AG, **IN** the last intron
        // rule 
        //      Variant in donor GT or acceptor AG, ***NOT IN***  the last intron, ....
        // then why holmes only check "NOT IN" last?
        if(trans.strand == '+'){
            for( size_t i{}; i<trans.exon.size()-1; ++i )
                if( trans.exon[i].first - pos <= 2 )
                    return fa.check( chr, trans.exon[i].first - 2, trans.exon[i].first - 1, "AG" );
        }else{
            for( size_t i{1}; i<trans.exon.size(); ++i )
                if( pos - trans.exon[i].second <= 2 ) 
                    return fa.check( chr, trans.exon[i].second + 1, trans.exon[i].second + 2, "CT" );
        }
        return false;
    }

    /**
     * @brief Checks if a variant affects the highly conserved donor GT sequence in a transcript.
     * 
     * This function checks if a variant at a given position in a transcript affects the highly
     * conserved donor GT sequence, which is located at the 5' end of introns and is essential for
     * proper splicing of pre-mRNA. The function only checks variants that are within 2 bases of
     * the donor GT sequence.
     * 
     * @param gene_id The Ensembl gene ID of the gene containing the transcript.
     * @param trans_id The Ensembl transcript ID of the transcript.
     * @param chr The chromosome where the variant is located.
     * @param pos The 0-based position of the variant in the transcript.
     * @param fa A reference to a `Fasta` object that provides access to the reference genome.
     * @return true if the variant affects the donor GT sequence, false otherwise.
     */
    bool check_donor(
        const std::string& gene_id, const std::string& trans_id,
        const std::string& chr, const size_t& pos, DataBaseFasta& fa, Variant& variant)
    {
        auto trans_opt = find_transcript(gene_id, trans_id);
        if(!trans_opt){
            // Maybe the hgvs is in RefSeq format
            // Use quick workaround by directly finding '+1G' or '+2T' in HGVSc
            return 
                (variant.hgvsc.find("+1G") != std::string::npos) or
                (variant.hgvsc.find("+2T") != std::string::npos);
        }
        decltype(auto) trans = trans_opt.value().get();
        
        // Check if the variant position is in the highly conserved donor GT 
        // |******exon****|GT-------intron----- ...
        //                 ^^ (these two position)
        // TODO: But!!!! Same as acceptor
        // rule EV0198:
        //      Variant in donor GT or acceptor AG, **IN** the last intron
        // rule 
        //      Variant in donor GT or acceptor AG, ***NOT IN***  the last intron, ....
        // then why holmes only check "NOT IN" last?
        if(trans.strand == '+'){
            for( size_t i{}; i<trans.exon.size()-1; ++i )
                if( pos - trans.exon[i].second <= 2 )
                    return fa.check( chr, trans.exon[i].second + 1, trans.exon[i].second + 2, "GT" );
        }else{
            for( size_t i{1}; i<trans.exon.size(); ++i )
                if( trans.exon[i].first - pos <= 2 )
                    return fa.check( chr, trans.exon[i].first - 2, trans.exon[i].first - 1, "AC" );
        }
        return false;
    }

    /**
     * @brief Check if the variant interrupt the last nucleotide G of the exon.
     *  Example: |*****exon******G|-------intron---- ....
     *                           ^ this position has variant
     * @param gene_id The Ensembl gene ID of the gene containing the transcript.
     * @param trans_id The Ensembl transcript ID of the transcript.
     * @param chr The chromosome where the variant is located.
     * @param pos The 0-based position of the variant in the transcript.
     * @param fa A reference to reference genome
     * @return true if the variant interrupts the last nucleotide G of an exon, false otherwise. 
     */
    bool check_splice_exon( const std::string& gene_id, const std::string& trans_id, const std::string& chr, const size_t& pos, DataBaseFasta& fa )
    {
        auto trans_opt = find_transcript(gene_id, trans_id);
        if(!trans_opt)
            return false;
        decltype(auto) trans = trans_opt.value().get();
    
        for(auto [pri5, pri3]: trans.exon){ // 5' pos and 3' pos of exons
            auto last_nucl_pos = (trans.strand == '+' ? pri3 : pri5);
            if(pos == last_nucl_pos)
                return fa.check_base( chr, last_nucl_pos,
                    (trans.strand == '+' ? 'G' : 'C'));
        }
        return false;
    }

    bool check_splice_intron(
        const std::string& gene_id, const std::string& trans_id,
        const std::string& chr, const size_t& pos, DataBaseFasta& fa, Variant& variant)
    {
        auto trans_opt = find_transcript(gene_id, trans_id);
        if(!trans_opt){
            // Maybe the hgvs is in RefSeq format
            // Use quick workaround by directly finding blablabla in HGVSc
            auto matching = [&hgvsc = variant.hgvsc](const auto& pattern){
                return hgvsc.find(pattern) != std::string::npos;
            };
            return matching("+3A") or matching("+3G") or matching("+4A") or matching("+5G");
        }
        decltype(auto) trans = trans_opt.value().get();
    
        // check if the variant:
        //  (1) located at the +3, +4 or +5 position of the intron
        //  (2) within the pathogenic range (i.e. <9 alleles in ExAC)
        //  (3) the reference nucleotide is:
        //      - Donor [+3 A/G] "OR" [+4 A] "OR" [+5 G] -> 2P (EV0184)
        // |*****exon*******|--AAG------intron---- ....
        //                     GAG
        //                     ^^^ variant involves in this range;
        for( auto [pri5, pri3]: trans.exon ){
            auto diff = (trans.strand == '+' ? pos - pri3 : pri5 - pos);
            if(diff >= 3 and diff <= 5)
                return trans.strand == '+' ?
                    fa.check_base(chr, pri3 + 3, 'A') || fa.check_base(chr, pri3 + 3, 'G') ||
                    fa.check_base(chr, pri3 + 4, 'A') ||
                    fa.check_base(chr, pri3 + 5, 'G')
                    :
                    fa.check_base(chr, pri5 - 3, 'T') || fa.check_base(chr, pri5 - 3, 'C') ||
                    fa.check_base(chr, pri5 - 4, 'T') ||
                    fa.check_base(chr, pri5 - 5, 'C');

                    // fa.check( chr, pri3 + 3, pri3 + 5, "AAG" ) || fa.check( chr, pri3 + 3, pri3 + 5, "GAG" ) : // + strand
                    // fa.check( chr, pri5 - 5, pri5 - 3, "CTT" ) || fa.check( chr, pri5 - 5, pri5 - 3, "CTC" );  // - strand

            // TODO: this is the original implement, is AGG/GGG correct?
            // is 'OR' (+3A/G or +4A or +5G) or 'AND' (AAG/GAG)
            // if( pos - i.second >=3 && pos - i.second <= 5 )  return fa.check( chr, i.second +3, i.second +5, "AGG" ) || fa.check( chr, i.second +3, i.second +5, "GGG" );
        }
        return false;
    }

    [[gnu::deprecated]] bool is_in_last_intron(const std::string& gene_id, const std::string& trans_id, const std::string& chr, const size_t& pos){
        auto trans_opt = find_transcript(gene_id, trans_id);
        if(!trans_opt)
            return false;
        decltype(auto) trans = trans_opt.value().get();

        if(trans.exon.size() < 2){
            return false;
        }

        using ExonType = std::pair<int, int>;
        auto last_idx = trans.exon.size() - 1;
        auto last2_exons = (trans.strand == '+') ?
            std::pair<ExonType, ExonType>{trans.exon[last_idx-1], trans.exon[last_idx]} :
            std::pair<ExonType, ExonType>{trans.exon[0]         , trans.exon[1]       };

        return pos > last2_exons.first.second and pos < last2_exons.second.first;
    }
/*
    int check_cnv( const std::string& gene_id, const std::string& trans_id, const size_t& start, const size_t& end )
    {
        auto trans = gene_map.find( gene_id )->second.trans_map.find( trans_id )->second;
        if( trans.exon[0].first > trans.exon[1].first ) std::reverse( trans.exon.begin(), trans.exon.end() );
        double s(0.5), e(0.5);
        size_t max = trans.exon.size();
        for( size_t i = 0; i<max; ++i )
        {
            if( start > trans.exon[i].first && start < trans.exon[i].second )  s = i+1;
            if( i != max && start > trans.exon[i].second && start < trans.exon[i+1].first ) s = i + 1 + 0.5;
            if( end > trans.exon[i].first && end < trans.exon[i].second ) e = i + 1;
            if( i != max && end > trans.exon[i].second && end < trans.exon[i+1].first )  e = i + 1 + 0.5;
            if( start > trans.exon[max-1] ) s = max + 0.5;
            if( end > trans.exon[max-1] ) e = max + 0.5;
        }
        if( s<=1 && e>=max )  return 1;
        else if( s<=1 )  return 2;
        else if( e>=max ) return 3;
    }
*/    
    int check_pos( const std::string& gene_id, const std::string& trans_id, const size_t& pos )
    {
        // TODO: This funtion is never called, and the purpose is unclear
        // But I will still modify the code
        auto trans_opt = find_transcript(gene_id, trans_id);
        if(!trans_opt)
            return false;
        decltype(auto) trans = trans_opt.value().get();
        int num = 0;
        for( auto [pri5, pri3] : trans.exon ){
            num = (trans.strand == '+') ? pos - pri3 : pri5 - pos;
            if( num <= 5 and num > -1 )
                return num;
        }
        return -1;
    }

    auto location_at(const std::string& gene_id, const std::string& trans_id, const size_t& pos){
        auto trans_opt = find_transcript(gene_id, trans_id);
        if(!trans_opt)
            return std::make_tuple(' ', ' ', std::pair<int, int>{-1, -1});
        decltype(auto) trans = trans_opt.value().get();

        char feature = 'E'; // E for exon, I for Intron
        char strand = trans.strand;
        auto& exons = trans.exon;
        std::pair<int, int> feature_intv = {-1, -1};
        auto within_interval = [pos](const std::pair<int, int>& intv){
            return pos >= intv.first and pos <= intv.second;
        };

        if(within_interval(exons[0])){
            feature_intv = exons[0];
        }else{
            for(int i = 1; i < exons.size(); ++ i){
                // intron type
                auto intron = std::pair<int, int>(exons[i-1].second+1, exons[i].first-1);
                if(within_interval(intron)){
                    feature = 'I';
                    feature_intv = intron;
                    break;
                }else if(within_interval(exons[i])){
                    feature_intv = exons[i];
                    break;
                }
            }
        }
        return std::make_tuple(feature, strand, feature_intv);
    }

};

}
