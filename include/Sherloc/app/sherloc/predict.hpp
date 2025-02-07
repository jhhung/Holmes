#pragma once 

#include <vector>
#include <string>
#include <Sherloc/sherloc_member.hpp>
#include <Sherloc/Patient/patient.hpp>
#include <Sherloc/DB/dbset.hpp>
#include <spdlog/spdlog.h>
#include <cmath>

namespace Sherloc::app::sherloc {

class Predict
{
  public:
    enum CompPred{
        pd0,
        pd1,
        pd2,
        pd3,
        pd4,
        pd5,
        pd6,
        pd7
    };

    /**
     * @brief General protein predictors, currently: SIFT, Polyphen, REVEL, CADD
     * 
     * @param variant 
     * @param para 
     */
    void general_algorithm(
        Variant & variant,
        const SherlocParameter& para 
    ) {
        static constexpr double revel_threshold = 0.5;
        static constexpr double cadd_threshold = 25.;

        auto revel_score = variant.revel_score;
        if(std::isnan(revel_score))
            revel_score = 0.; // TODO: this can be discussed

        auto cadd_score = variant.cadd_phred_score;
        auto has_cadd = !std::isnan(cadd_score);
        if(variant.sift and variant.poly and revel_score >= revel_threshold){ // all deleterious
            variant.add_rule( 122 );
            variant.add_tag(HOLMES_MAKE_TAG(pd0));
            if(has_cadd and cadd_score >= cadd_threshold){
                variant.add_rule( 217 );
            }
        } else if(!variant.sift and !variant.poly and revel_score < revel_threshold){ // all benign
            variant.add_rule( 126 );
            variant.add_tag(HOLMES_MAKE_TAG(pd2));
            if(has_cadd and cadd_score < cadd_threshold){
                variant.add_rule( 218 );
            }
        } else { // conflict
            variant.add_rule( 109 );
            variant.add_tag(HOLMES_MAKE_TAG(pd1));
        }
    
    }

    /**
     * @brief Gene-specific predictors
     * 
     * @param variant 
     * @param para 
     */
    void gene_specific_algorithm(
        Variant & variant,
        const SherlocParameter& para 
    ) {
        // TODO: find some gene-specific tools
        // For example: PMID 19043619, 23621914, 18383312, 21310275, 22290698
        // need to include rule: EV0121, EV0127
    }

    /**
     * @brief Mammalian conservation alleles
     * 
     * @param variant 
     * @param para 
     */
    void mammalian_convervation(
        Variant & variant,
        const SherlocParameter& para 
    ) {
        // TODO: variant alleles that are conserved in mammalian species
        // need to include rule: EV0125
    }

    void go_protein( SherlocMember& sher_mem
            , const size_t& index 
            , const SherlocParameter& para 
            )
    {
        auto& variant = sher_mem.variants[index];
        general_algorithm(variant, para);
        gene_specific_algorithm(variant, para);
        mammalian_convervation(variant, para);
    }

    void go_splice_loss( SherlocMember& sher_mem
            , const size_t& index 
            , const SherlocParameter& para 
            )
    {
        if( sher_mem.variants[index].mes_score < 0 ){
            sher_mem.variants[index].add_rule( 187 );
            sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(pd3));
        }else{
            sher_mem.variants[index].add_rule( 191 );
            sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(pd7));
        }
    }

    void go_splice_gain( SherlocMember& sher_mem
            , const size_t& index 
            , const SherlocParameter& para 
            )
    {
        if( sher_mem.variants[index].mes_score >= 0.1 ){
            sher_mem.variants[index].add_rule( 202 );
            sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(pd6));
        }else if( sher_mem.variants[index].mes_score >= 0.03 ){
            sher_mem.variants[index].add_rule( 203 );
            sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(pd5));
        }else{
            sher_mem.variants[index].add_rule( 201 );
            sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(pd4));
        }
    }

    void go_splice(
        SherlocMember& sher_mem,
        size_t index,
        const SherlocParameter& para 
    ) {
        // mes is (MaxEntScan_diff < 0), and MaxEntScan_diff is (SWA_REF_COMP - SWA_ALT)
        // so when mes == true: score increased, mes == false: score decreased
        if( sher_mem.variants[index].mes == false )
            go_splice_loss( sher_mem, index, para );
        else
            go_splice_gain( sher_mem, index, para );
    }

    void go_multi_factorial(
        SherlocMember& sher_mem,
        size_t index,
        const SherlocParameter& para 
    ) {
        // TODO: Multi-factorial algorithm
        // need to include rule: EV0119, EV0118
    }

    void run( Patient::Patient& ze)
    {
        decltype(auto) para = SherlocParameter::get_paras();
        for( auto& sher_mem: ze.sher_mems )
        {
            for( size_t j{}; j<sher_mem.variants.size(); ++j )
            {
                bool is_missense = std::ranges::find(sher_mem.variants[j].type, "missense_variant") != sher_mem.variants[j].type.end();
                if(is_missense){
                    // Protein effect (missense changes only)
                    go_protein( sher_mem, j, para );
                }

                // TODO: currently we only use MES, so when MES is not empty, we go to splice rule
                // But there might be better condition
                bool is_splice = !sher_mem.variants[j].mes_empty;
                if(is_splice){
                    // Splice effect
                    go_splice(sher_mem, j, para);
                }

                go_multi_factorial(sher_mem, j, para);
            }
        }
    }
};

}
