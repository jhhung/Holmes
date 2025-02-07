#pragma once

#include <Sherloc/Patient/patient.hpp>
#include <Sherloc/disease_database.hpp>
#include <Sherloc/sherloc_member.hpp>
#include <Sherloc/app/sherloc/sherloc_parameter.hpp>
#include <Sherloc/DB/dbset.hpp>
#include <Sherloc/Patient/other_patient.hpp>

namespace Sherloc::app::sherloc {

class Clinical
{
  public:
    enum ClinicalSegregation{
        cs0,
        cs1,
        cs2,
        cs3
    };
    enum ClinicalTree1{
        ca0,
        ca1,
        ca2,
        ca3,
        ca4,
        ca5,
        ca6,
        ca7,
        ca8,
        ca9,
        ca10
    };
    enum ClinicalTree2{
        cb0,
        cb1,
        cb2,
        cb3,
        cb4,
        cb5,
        cb6,
        cb7,
        cb8,
        cb9,
        cb10,
        cb11,
        cb12,
        cb13,
        cb14
    };
    enum ClinicalTree3{
        cc0,
        cc1,
        cc2,
        cc3,
        cc4,
        cc5,
        cc6,
        cc7
    };

    void go_tree3( 
             const SherlocParameter& para
           , Patient::Patient& ze
         )
    {
        for( auto& sher_mem: ze.sher_mems )
        {
            auto genotype = sher_mem.genotype;

            // get the decision evidence
            // FIXME: Is my understanding right?
            auto is_AD_or_Xlinked_or_Ylinked =
                sher_mem.is_autosomal_dominant() or
                (sher_mem.is_x_linked() and ze.sex) or
                sher_mem.is_y_linked();

            auto is_early_onset = sher_mem.onset;

            auto is_genotype_homozygous = 
                (genotype[0] != -1 and genotype[1] != -1) and
                (genotype[0] == genotype[1]);

            
            if(is_AD_or_Xlinked_or_Ylinked){
                if(is_genotype_homozygous){
                    if(is_early_onset){
                        sher_mem.add_rule( 130 );
                        sher_mem.add_tag(HOLMES_MAKE_TAG(cc1));
                    }else{
                        sher_mem.add_rule( 84 );
                        sher_mem.add_tag(HOLMES_MAKE_TAG(cc0));
                    }
                }
                else{
                    if(is_early_onset){
                        sher_mem.add_rule( 134 );
                        sher_mem.add_tag(HOLMES_MAKE_TAG(cc2));
                    }else{
                        sher_mem.add_rule( 53 );
                        sher_mem.add_tag(HOLMES_MAKE_TAG(cc3));
                    }
                }
            }
            else{
                if(is_genotype_homozygous){
                    if(is_early_onset){
                        sher_mem.add_rule( 129 );
                        sher_mem.add_tag(HOLMES_MAKE_TAG(cc6));
                    }else{
                        sher_mem.add_rule( 84 ); 
                        sher_mem.add_tag(HOLMES_MAKE_TAG(cc7));
                    }
                }
                // FIXME: these two look like rules that need compound hetero
                // else if( sher_mem.next != 0 && (*(&sher_mem + sher_mem.next)).consequence == true ){
                //     if(is_early_onset)
                //         sher_mem.add_rule( 140 );
                //     else  
                //         sher_mem.add_rule( 141 );
                // }
                else continue;
            }
            sher_mem.clinical_rule = true;
        }
    }

    void go_greater(
             const SherlocParameter& para
           , Patient::Patient& ze
           , Patient::OtherPatient& op   
        )
    {
        for(auto idx = 0; idx < ze.sher_mems.size(); ++idx)
        {
            auto& sher_mem = ze.sher_mems[idx];
            if (sher_mem.af_above_somewhat_high()) {
                // According to Supplementary file CASE TREE #1
                // " Variant Frequency: Somewhat high", "Segregation analysis only"
                go_family(para, sher_mem, ze);
                continue;
            }
            int op_num = op.get_observation( sher_mem )+1;

            auto denovo_status = ze.check_allele_denovo(idx);
            if (denovo_status == Attr::Allele::IsDeNovo) {
                sher_mem.add_rule( 205 );
            } else if (denovo_status == Attr::Allele::Unknown) {
                if (ze.is_denovo) {
                    sher_mem.add_rule( 205 );
                } else {
                    // There are some constraints:
                    // [v] "... If the variant is common (above somewhat high MAF), 
                    // [x]  the gene with the de novo variant is not well-established to cause disease, TODO: 
                    // [v]  the disease is low penetrance, <- This branch is "go_greater", so it's already assumed to be high penetrance
                    // [v]  and/or parents are affected with disease, 
                    //  do not apply this criteria"
                    if(!ze.is_dad_sick and !ze.is_mom_sick)
                        sher_mem.add_rule( 206 );
                }
            }

            for(int var_idx = 0; var_idx < sher_mem.variants.size(); ++var_idx){
                auto& var = sher_mem.variants[var_idx];
                auto inhe_patt = sher_mem.inheritance_patterns[var_idx];
                
                if (inhe_patt == 'D' or
                    (inhe_patt == 'X' and ze.sex) or
                    inhe_patt == 'Y') {
                    if (sher_mem.gt_hetero() or sher_mem.gt_unknown()){
                        for( int i{}; i<op_num; ++i )
                            var.add_rule( 169 ); 
                    }
                } else if (
                    inhe_patt == 'R' or
                    (inhe_patt == 'X' and not ze.sex) or
                    inhe_patt == 'U'
                ) {
                    if(sher_mem.gt_homo()){
                        var.add_rule( 153 ); // TODO: there are alot of contraints
                        var.add_tag(HOLMES_MAKE_TAG(ca0));
                    }
                }
            }
        }
    }

    void case_report(Variant& var, int observation_patient_num){
        if( observation_patient_num < 2 )
            var.add_rule( 211 );
        else if( observation_patient_num < 3 )
            var.add_rule( 81 );
        else if( observation_patient_num < 4 )
            var.add_rule( 80 );
        else
            var.add_rule( 79 );
    }

    void go_less(
             const SherlocParameter& para
           , Patient::Patient& ze
           , const double& yield
           , Patient::OtherPatient& op   
        )
    {
        for(auto idx = 0; idx < ze.sher_mems.size(); ++idx)
        {
            auto& sher_mem = ze.sher_mems[idx];
            if (sher_mem.af_above_somewhat_high()) {
                // According to Supplementary file CASE TREE #1
                // " Variant Frequency: Somewhat high", "Segregation analysis only"
                go_family(para, sher_mem, ze);
                continue;
            }

            auto dad = ze.dad_vcf.find( sher_mem );
            auto mom = ze.mom_vcf.find( sher_mem );
            auto dad2 = dad;
            auto mom2 = mom;
            if( sher_mem.next != 0 )
            {
                dad2 = ze.dad_vcf.find( *(&sher_mem + sher_mem.next) );
                mom2 = ze.mom_vcf.find( *(&sher_mem + sher_mem.next) );
            }
            int op_num = op.get_observation( sher_mem )+1;

            auto denovo_status = ze.check_allele_denovo(idx);
            if (denovo_status == Attr::Allele::IsDeNovo) {
                sher_mem.add_rule( 205 );
            } else if (yield >= 15 and denovo_status == Attr::Allele::Unknown) { // not low penetration
                if (ze.is_denovo) {
                    sher_mem.add_rule( 205 );
                } else {
                    // There are some constraints:
                    // [v] "... If the variant is common (above somewhat high MAF), 
                    // [x]  the gene with the de novo variant is not well-established to cause disease, TODO: 
                    // [?]  the disease is low penetrance, <- TODO: How low? we use < 15 for now
                    // [v]  and/or parents are affected with disease, 
                    //  do not apply this criteria"
                    if(!ze.is_dad_sick and !ze.is_mom_sick)
                        sher_mem.add_rule( 206 );
                }
            }
            
            // for each variant
            for(int var_idx = 0; var_idx < sher_mem.variants.size(); ++var_idx){
                auto& var = sher_mem.variants[var_idx];
                auto inhe_patt = sher_mem.inheritance_patterns[var_idx];
                
                if(
                    (inhe_patt == 'D' or (inhe_patt == 'X' and ze.sex == true))
                ){
                    if (sher_mem.gt_hetero() or sher_mem.gt_unknown()) { // unknown is assume to be hetero)
                        case_report(var, op_num);
                    }
                } else if (inhe_patt == 'R' or inhe_patt == 'U') { // assume unknown inheritance pattern as recessive
                    if (sher_mem.gt_hetero() or sher_mem.gt_unknown()){
                        // case1: 1 variant or 2 variants in cis
                        //      '2 variants in cis' is not handled here
                        var.add_rule( 107 );
                    } else {
                        // case2: 2 variants, phase unknown
                        //      not handled here

                        // case3: 2 in trans or homoygous
                        case_report(var, op_num);
                    }
                }
            }
        }
    }

    void go_tree1(
             const SherlocParameter& para
           , Patient::Patient& ze
           , const double& yield
           , Patient::OtherPatient& op   
          )
    {
        if( yield > 75 )
            go_greater( para, ze, op );
        else
            go_less( para, ze, yield, op );
    }

    void go_tree2(){
        // TODO: some benign rules, need "explanatory variant" in same / another gene
        //   and phase information, is hard to implement in compound_heterozygous_candidate, too.
        // EV0132, EV0133, EV0060, EV0061
    }

    void go_family(
             const SherlocParameter& para
           , SherlocMember& sher_mem
           , Patient::Patient& patient
        //    , std::vector< Patient::Patient >& patient 
          )
    {
        // TODO: family segregation should be rewritten when we got some family data for validation or testing

        // for( auto& ze: patient )
        //     for( auto& sher_mem: ze.patient )
        //     {
        //         int healthy{}, sick{};
        //         for( auto& ze2: patient ) // FIXME: so weird... what the purpose of nested for loop through other patients? this is "Segregation", how can it be "Segregation" if you count other patient's families?
        //         {
        //             sick += ze2.sick_family.get_observation( sher_mem );
        //             healthy += ( ze2.healthy_family.op_vec.size() - ze2.healthy_family.get_observation( sher_mem ) );
        //         }
        //         if( sher_mem.is_autosomal_dominant() )
        //         {
        //             if( sick < 3 )  continue;
        //             if( sick < 6  ) sher_mem.add_rule( 51 );
        //             if( sick < 10 ) sher_mem.add_rule( 50 );
        //             if( sick >=10 ) sher_mem.add_rule( 49 );
        //             sher_mem.clinical_rule = true;
        //         }
        //         else
        //         {
        //             if( sick < 1 )  continue;
        //             if( sick < 2 && healthy < 3 ) continue;
        //             if( sick < 3 && healthy < 4 )  sher_mem.add_rule( 51 );
        //             if( sick < 3 || sick < 5 && healthy < 7 )  sher_mem.add_rule( 50 );
        //             else sher_mem.add_rule( 49 );
        //             sher_mem.clinical_rule = true;
        //         }
        //     }
        // }
    }

    void run(
           Patient::Patient& patient 
           , const Disease& disease
           , Patient::OtherPatient& op   
          )
    {
        decltype(auto) para = SherlocParameter::get_paras();
        if (patient.is_sick) { // Observations of variants in affected individuals
            // TODO: "Has an alternate cause of disease"
            // Have not implemented yet. Currently we only go tree 1
            bool alternative_cause = false;

            if (!alternative_cause) {
                go_tree1( para, patient, disease.yield, op );
            } else {
                go_tree2();
            }
        } else { // Observations in unaffected individuals
            go_tree3( para, patient );
        }
    }
};

}
