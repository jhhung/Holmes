#pragma once

#include <vector>
#include <set>
#include <string>
#include <Sherloc/DB/gtf.hpp>
#include <Sherloc/Patient/patient.hpp>
#include <Sherloc/app/sherloc/sherloc_parameter.hpp>
#include <Sherloc/DB/uniprot.hpp>
#include <Sherloc/app/sherloc/sherloc_consequence.hpp>
#include <Sherloc/DB/fasta.hpp>
#include <Sherloc/DB/dbset.hpp>

namespace Sherloc::app::sherloc {

class VariantRule {
  public:
    // TODO: maybe construct a SE (sequence ontology) graph representing all variant types
    // reference: https://asia.ensembl.org/info/genome/variation/prediction/predicted_data.html#consequences
    std::set<std::string> null_types;
    std::set<std::string> splice_types;
    std::set<std::string> missense_types;
    std::set<std::string> initiator_types;
    std::set<std::string> silent_types;
    std::set<std::string> intronic_types;
    std::set<std::string> inframe_types;
    std::set<std::string> noncoding_types;

    VariantRule():
        null_types{
            "stop_gained",
            "frameshift_variant"
        },
        splice_types{
            "splice_acceptor_variant",
            "splice_donor_variant",
            "splice_region_variant",
            "splice_donor_region_variant",
            "splice_donor_5th_base_variant"
        },
        missense_types{
            "missense_variant"
        },
        initiator_types{
            "start_lost"
        },
        silent_types{
            "synonymous_variant"
        },
        intronic_types{
            "intron_variant"
        },
        inframe_types{
            "inframe_deletion",
            "inframe_insertion"
        },
        noncoding_types{
            "5_prime_UTR_variant",
            "3_prime_UTR_variant",
            "intergenic_variant",
            "upstream_gene_variant",
            "downstream_gene_variant",
            "non_coding_transcript_exon_variant",
            "non_coding_transcript_variant"
        }
    {}

    enum VarTypeInframe{
        vi0
    };
    enum VarTypeInitiator{
        vin0,
        vin1,
        vin2
    };
    enum VarTypeIntronic{
        vt0,
        vt1
    };
    enum VarTypeMissense{
        vm0,
        vm1,
        vm2,
        vm3,
        vm4
    };
    enum VarTypeNoncoding{
        vo0
    };
    enum VarTypeNull{
        vn0,
        vn1,
        vn2,
        vn3,
        vn4,
        vn5,
        vn6,
        vn7
    };
    enum VarTypeSilent{
        vsi0,
        vsi1
    };
    enum VarTypeSplice{
        vs0,
        vs1,
        vs2,
        vs3,
        vs4,
        vs5,
        vs6,
        vs7
    };

    inline bool loss_of_function(
        SherlocMember& sher_mem,
        int variant_index
    ){
        static auto is_rule_which_lof = [](auto rule){
            // FIXME: I just follow what jeff written, but why a variant is lof with these rules?? 
            return rule == 50 or rule == 49 or rule == 80 or rule == 79;
        };
        bool lof = false;
        if( sher_mem.variants[variant_index].pli >= 0.9 ){
            lof = true;
        }else{
            lof = std::ranges::any_of(sher_mem.group, is_rule_which_lof);
        }
        return lof;
    }

    inline bool is_snv(const DB::Clinvar& var){
        return 
            var.ref != "-" and 
            var.alt != "-" and 
            var.ref.size() == 1 and 
            var.alt.size() == 1;
    };

    void go_null( SherlocMember& sher_mem
                , const int& index
                , DB::DataBaseGTF& gtf
                , const SherlocParameter& para
                , const DB::DBSet& db
                , SherlocConsequence& sher_conseq  
            ) 
    {
        bool lof = loss_of_function(sher_mem, index);

        decltype(auto) variant = sher_mem.variants[index];

        bool nmd = std::ranges::any_of(variant.type,
            [](auto& var_type){
                return var_type == "NMD_transcript_variant";
            });

        if(lof){ // Loss-of-function established
            if(nmd and !variant.might_escape_nmd){
                variant.add_rule( 16 );
                variant.add_tag(HOLMES_MAKE_TAG(vn0));
            }else if( sher_conseq.get_null( variant.trans, sher_mem.pos ) == true ){
                variant.add_rule( 175 );
                variant.add_rule( 19 );
                variant.add_tag(HOLMES_MAKE_TAG(vn4));
            }else if( sher_conseq.get_miss( variant.trans, sher_mem.pos ) == true ){
                variant.add_rule( 19 );
                variant.add_rule( 44 );
                variant.add_tag(HOLMES_MAKE_TAG(vn3));
            }else{
                variant.add_rule( 19 );
                variant.add_tag(HOLMES_MAKE_TAG(vn2));
            }
        }else{ // Loss-of-function NOT established (can include gain-of-function)
            if(nmd and !variant.might_escape_nmd){
                variant.add_rule( 183 );
                variant.add_tag(HOLMES_MAKE_TAG(vn1));
            }else if( sher_conseq.get_null( variant.trans, sher_mem.pos ) == true ){
                variant.add_rule( 175 );
                variant.add_rule( 194 );
                variant.add_tag(HOLMES_MAKE_TAG(vn5));
            }else{
                variant.add_rule( 194 );
                variant.add_tag(HOLMES_MAKE_TAG(vn6));
            }
        }
    }

    void go_splice( SherlocMember& sher_mem
                  , const int& index
                  , DB::DataBaseGTF& gtf
                  , const SherlocParameter& para
                  , DB::DataBaseFasta& fa
                ) 
    {
        bool lof = loss_of_function(sher_mem, index);
           
        bool is_frameshift(false), is_donor(false), is_acceptor(false);

        decltype(auto) variant = sher_mem.variants[index];
        for( auto& i: variant.type )
        {
            // All var types that go to this tree:
            // "splice_acceptor_variant", "splice_donor_variant",
            // "splice_region_variant", "splice_donor_region_variant",
            // "splice_donor_5th_base_variant"

            // TODO: correct me if I'm wrong: 
            // I don't think the image of sherloc supplementary `gim201737x4.pdf`
            // splice variant: "... affected exon disrupts reading frame" means this variant has frame_shift type
            // it just means that the variant satisfy one of:
            //   1) +GT/-AG and not in last intron
            //   2) Last nucleotide of exon (G only)
            //   3) Donor +3 A/G, +4A, +5G (canonical +GT junction only)
            // if( i == "frameshift_variant" ) is_frameshift = true; 
            if( i == "splice_acceptor_variant" )  is_acceptor = true;
            if( i == "splice_donor_variant" ) is_donor = true;
        }

        auto donor_GT_or_acceptor_AG = 
            (is_donor and gtf.check_donor( variant.gene, variant.trans, sher_mem.chr, sher_mem.pos, fa, variant))
            or
            (is_acceptor and gtf.check_acceptor( variant.gene, variant.trans, sher_mem.chr, sher_mem.pos, fa, variant));

        // DEPRECATED, use information from VEP '--numbers' option annotated INTRON is easier
        // auto in_last_intron = gtf.is_in_last_intron(variant.gene, variant.trans, sher_mem.chr, sher_mem.pos);

        // New way that use INTRON information
        auto in_last_intron = (
            variant.sub_feat.has and 
            !variant.sub_feat.is_exon and 
            variant.sub_feat.idx == variant.sub_feat.total
        );
        SPDLOG_DEBUG("in_last_intron: {}, {}/{}", in_last_intron,
            variant.sub_feat.has_value() ? variant.sub_feat->idx : 0,
            variant.sub_feat.has_value() ? variant.sub_feat->total : 0);


        if(lof){
            if( gtf.check_splice_exon( variant.gene, variant.trans, sher_mem.chr, sher_mem.pos, fa ) ){
                variant.add_rule( 196 );
                variant.add_tag(HOLMES_MAKE_TAG(vs3));
            }
            if( gtf.check_splice_intron( variant.gene, variant.trans, sher_mem.chr, sher_mem.pos, fa, variant) ){
                variant.add_rule( 184 );
                variant.add_tag(HOLMES_MAKE_TAG(vs4));
            }
            if(donor_GT_or_acceptor_AG){
                if(in_last_intron){
                    variant.add_rule( 198 );
                    variant.add_tag(HOLMES_MAKE_TAG(vs1));
                }else{
                    variant.add_rule( 17 );
                    variant.add_tag(HOLMES_MAKE_TAG(vs2));
                }
            }
        }else{
            if(donor_GT_or_acceptor_AG and !in_last_intron){
                variant.add_rule( 181 );
                variant.add_tag(HOLMES_MAKE_TAG(vs5));
            }
        }
    }

    void go_missense( SherlocMember& sher_mem
                  , const int& index
                  , const SherlocParameter& para
                  , DB::DataBaseUniprot& uniprot
                  , DB::DataBaseGTF& gtf
                  , DB::DataBaseClinvar& clinvar
                ) 
    {
        auto& variant = sher_mem.variants[index];
        // New impl

        // EV0018  P       4.0     "Same AA change as a different pathogenic variant"
        // "Usage notes": 
        //    Use for a novel nucleotide change that results in the 
        //    [[exact same amino acid change]] as a previously documented pathogenic missense change. 
        //    Do not use if the mechanism of pathogenicity for the reported mutation is a splicing defect. 
        //                                                           (How to know??)      ^
        //                                                     TODO: Maybe from clinvar INFO MC entry, (e.g. MC=SO:0001575|splice_donor_variant)
        //    Must first determine that the co-occurring variant meets our [[criteria for Pathogenic or Likely Pathogenic]]. 
        //    To be pathogenic, this variant type should be novel or very rare; 
        //    additional points for EV0135 or EV0101 will not be added to final calculated score. 
        //    If the allele count for the variant is above the expected pathogenic range do not apply this critiera. 
        //    Instead, apply clinical findings and experimental studies evidence as appropriate.
        
        // EV0044  P       1.5     "Missense change at ""important"" amino acid residue"
        //    An ""important"" amino acid is inferred from previous reports of pathogenic 
        //    (should meet our criteria for LP or P) missense changes at this codon. 
        //    To use this critiera, the previously reported missense change must meet 
        //    Sherloc guidelines for Pathogenic or Likely Pathogenic classifications. 
        //    (If the mechanism of pathogenicity for the previously reported missense 
        //    change is known to be splicing, don't use this criteria.)
        bool use_rule44_if_no_rule18 = false;


        // New approach with DB_Clinvar
        auto reported_variants_in_same_aa = 
            clinvar.find_txp_by(
                variant.trans,
                {variant.aa_pos, variant.aa_pos},
                &DB::DataBaseClinvar::AA_proj
            );

        for(auto v : reported_variants_in_same_aa){
            if(v and is_snv(*v)){
                auto clinvar_var_codon = uniprot.get_codon( v->codon );
                if(v->consequence){
                    if(variant.codon == clinvar_var_codon){ // [[exact same amino acid change]]
                        variant.add_rule( 18 );
                        variant.add_tag(HOLMES_MAKE_TAG(vm0));
                        use_rule44_if_no_rule18 = false;
                        break;
                    }else{
                        SPDLOG_CRITICAL("Clinvar Allele: ref{}, alt{}", v->ref, v->alt);
                        SPDLOG_CRITICAL("Rule 44: clinvar id {} codon: {}, uniprot: {}", 
                            v->allele_id, clinvar_var_codon, v->codon);
                        use_rule44_if_no_rule18 = true; // different missense changes at this codon
                    }
                }

                if(v->benign){ // FIXME: maybe need to check the AF, should above high
                    if(
                        // FIXME: should also be a "different nucleotide substitution"
                        // current implementation of DBClinvar is kinda poor
                        // can't get the exact position when query by CDS / AA pos
                        variant.codon == clinvar_var_codon 
                    ){ // [[exact same amino acid change]]
                        SPDLOG_DEBUG("var codon '{}' vs clinvar codon '{}'", variant.codon, clinvar_var_codon);
                        variant.add_rule( 112 );
                        variant.add_tag(HOLMES_MAKE_TAG(vm1));
                        break;
                    }
                }
            }
        }
        if(use_rule44_if_no_rule18){
            variant.add_rule( 44 );
            variant.add_tag(HOLMES_MAKE_TAG(vm3));
        }

        // EV0212  P       4.0     "Missense change at ""critical"" aa residue"
        // "Usage notes": 
        //    Critical residues are residues that are intolerant of any change;
        //    any alteration from reference completely disrupts protein structure and function.
        //    [[A rigorous method for establishing a residue as “critical” is not provided here and caution should be used.]]
        //    This criteria is generally only applied when substantial clinical evidence 
        //    about related residues in the same gene is available, and knowledge of 
        //    protein structure indicates which disparate residues are united by known essential physical interactions.
        // TODO: need to find some sources of "Critical residues"
        ///////// possible implementation:
        // if(critical_residues.contains(variant)){
        //     ...
        // }

        // EV0172  P       1.0     "Mutation rich region"
        // "Usage notes": 
        //    Used when multiple pathogenic missense changes have been reported 
        //    surrounding this missense change. The missense change in question should
        //    involve a highly conserved residue and the reported pathogenic missense 
        //    changes should involve similarly conserved residues. Without experimental support,
        //    an absolute minimum of reported mutations is not appropriate, but a good
        //    example of a mutation rich region is the Mitofusin HR2 domain of MFN2, 
        //    with 19 reported missense mutations between AAs 586 and 755 (see HGMD). 
        //    Also, experimental evidence supporting the functional importance of a domain
        //    (as is the case with HR2 of MFN2) is beneficial, but not sufficient in the 
        //    absence of reported pathogenic mutations. The missense changes should be classified
        //    as P or LP to be counted, although for very rare, highly penetrant diseases, 
        //    absent from ExAC and reported in at least 1 patient with expected disease is sufficient.
        // Need source for "Mutation rich region", Currently Clinvar, TODO: Uniprot for functional domain
        auto reported_variants_in_same_region = 
            clinvar.find_txp_by(
                variant.trans,
                {
                    (para.mut_rich_config.windows > variant.cds_pos ?
                        0 : variant.cds_pos - para.mut_rich_config.windows), // max(0, region start)
                    variant.cds_pos + para.mut_rich_config.windows
                },
                &DB::DataBaseClinvar::CDS_proj
            );

        int patho_count = 0, benign_count = 0;
        for(auto v : reported_variants_in_same_region){
            if(v and is_snv(*v)){ // only use for missense change
                if(v->consequence) ++patho_count;
                if(v->benign) ++benign_count;
            }
        }


        // TODO: current approach, could be changed if any better criteria exists
        //  1) greater than or equal to 4 pathogenic variants
        //  2) no benign variant
        if(
            patho_count >= para.mut_rich_config.pathogenic_variant_threshold and 
            benign_count == 0
        ) {
            variant.add_rule(172);
            // variant.add_tag(HOLMES_MAKE_TAG(vm5));
        }
    }

    void go_initiator( SherlocMember& sher_mem
                  , const int& index
                  , const SherlocParameter& para
                  , DB::DataBaseUniprot& uniprot
                ) 
    {
        bool lof = loss_of_function(sher_mem, index);

        if(!lof)
        {
            sher_mem.variants[index].add_rule( 182 );
            sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(vin0));
            return;
        }
        sher_mem.variants[index].add_rule( 114 );
        sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(vin1));
        // TODO: EV0018 can be added here
        // "... This can be added with EV0018: Same AA change as a pathogenic missense 
        // if a different variant within this start codon is reported as pathogenic. ..."
        // sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(vin2));
    }

    void apply_103_with_check(
        SherlocMember& sher_mem,
        const int& index,
        DB::DataBaseGTF &gtf
    ){
        // EV0103  B       2.0     "Silent and intronic changes outside of the consensus splice sites"
        // TODO: Use for a silent or intronic variant, except when it is located: 
        //  (1) at the -3 to -1 nucleotides within the intron,
        //  (2) at the +1 to +6 nucleotides within the intron,
        //  (3) at the first nucleotide of the exon, or
        //  (4) at the last 2 nucleotides of the exon.
        decltype(auto) variant = sher_mem.variants[index];
        auto within_interval = [pos = sher_mem.pos](const std::pair<int, int>& intv){
            return pos >= intv.first and pos <= intv.second;
        };

        auto [feature, strand, intv] = gtf.location_at(variant.gene, variant.trans, sher_mem.pos);
        if(feature != ' '){
            if(feature == 'I'){ // located at intron
                if((strand == '+') ?
                    within_interval({intv.second - 2, intv.second   }) :
                    within_interval({intv.first     , intv.first + 2})) // -3 to -1 nucleotides
                    return;

                if((strand == '+') ?
                    within_interval({intv.first     , intv.first + 5}) :
                    within_interval({intv.second - 5, intv.second   })) // +1 to +6 nucleotides
                    return;
            }else{
                if( (strand == '+' and sher_mem.pos == intv.first)
                    or
                    (strand == '-' and sher_mem.pos == intv.second)) // first nucleotide of the exon
                    return;

                if((strand == '+') ?
                    within_interval({intv.second - 1, intv.second   }) :
                    within_interval({intv.first     , intv.first + 1})) // last 2 nucleotides of the exon
                    return;
            }
        }else{
            // TODO: Can't find feature, might because of refseq cache
            // still need to find a way to approximate the location
            // currently still apply 103
        }
        sher_mem.variants[index].add_rule( 103 );
    }

    void go_silent(
        SherlocMember& sher_mem,
        const int& index,
        DB::DataBaseGTF &gtf
    ){
        apply_103_with_check(sher_mem, index, gtf);
        sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(vs1));
    }
    
    void go_intronic(
        SherlocMember& sher_mem, 
        const int& index, 
        DB::DataBaseGTF &gtf
    ){
        bool same = true;

        // TODO: check if the following homopolymer code is done correctly
        if( sher_mem.ref.size() > 1 )
            for( int i=1; i<sher_mem.ref.size(); ++i )
            {
                same = ( sher_mem.ref[i-1] == sher_mem.ref[i] );
                if( same == false ) break;
            }
        if( sher_mem.alt.size() > 1 && same == true )
            for( int i=1; i<sher_mem.ref.size(); ++i )
            {
                same = ( sher_mem.ref[i-1] == sher_mem.ref[i] );
                if( same == false ) break;
            }
        if( sher_mem.alt.size() + sher_mem.ref.size() > 2 && same == true ){
            sher_mem.variants[index].add_rule( 21 );
            sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(vt0));
        }else{
            apply_103_with_check(sher_mem, index, gtf);
            sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(vt1));
        }
    }

    void go_inframe( SherlocMember& sher_mem
                  , const int& index
                  , const SherlocParameter& para
                  , SherlocConsequence& sher_conseq
                )
    {
        sher_mem.variants[index].add_rule( 185 );
        sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(vi0));
        // EV0044  P       1.5     "Missense change at ""important"" amino acid residue"
        //    ... This criteria can also be used for small in-frame deletions, 
        //    single or multi exon deletions that encompass a known pathogenic missense, 
        //    and truncating mutations that are expected to escape NMD (EV0019). 
        //    To use this criteria with a truncating variant that may escape NMD,
        //          (1) the disease mechanism must be LoF-established and
        //          (2) the PTC created by the truncating variant must be upstream of the LP/P missense variant. 
        //    Must first determine that the co-occurring variant meets our criteria for Pathogenic or Likely Pathogenic.
        // TODO: this is complicated, should be carefully applied.
    }

    void go_non_coding( SherlocMember& sher_mem
                  , const int& index
                  , const SherlocParameter& para
                )
    {
        sher_mem.variants[index].add_rule( 168 );
        sher_mem.variants[index].add_tag(HOLMES_MAKE_TAG(vo0));
    } 

    void run( 
        SherlocMember& sher_mem
        , const SherlocParameter& para
        , DB::DBSet& db
        , SherlocConsequence& sher_conseq
    ) {
        // EV0139   P       1.5     "Disruption of essential nucleotide
        //    An ""essential"" nucleotide is inferred from previous reports of 
        //    pathogenic nucleotide changes at this same position. To use this critiera,
        //    the previously reported nucleotide change must meet Sherloc guidelines for 
        //    Pathogenic or Likely Pathogenic classifications. 
        //    For example, if this variant is a +5 G>C change , and a +5G>A change is 
        //    classified as Likely Pathogenic, this criteria would apply
        // So we give sher_mem a 139 rule if the variant is snv and the exact 
        // same pos is reported in ClinVar
        auto sher_mem_is_snv = (
            sher_mem.ref != "-" and 
            sher_mem.alt != "-" and
            sher_mem.ref.size() == 1 and
            sher_mem.alt.size() == 1
        );
        if (sher_mem_is_snv){
            const auto& clinvar = db.db_clinvar;
            auto it_chr = clinvar.chr2vec.find(Attr::ChrMap::chr2idx(sher_mem.chr));
            if (it_chr != clinvar.chr2vec.end()){
                auto& vec = it_chr->second;

                auto [s_it, e_it] = std::ranges::equal_range(
                    vec,
                    sher_mem.pos,
                    {},
                    &DB::DataBaseClinvar::PosClinvar::first
                );

                // check if any pathogenic variant at the same position in ClinVar is also SNV
                for(auto it = s_it; it != e_it; ++it){
                    auto& clinvar_variant = it->second;
                    if (is_snv(clinvar_variant) and clinvar_variant.consequence){
                        sher_mem.add_rule(139);
                        break;
                    }
                }
            }
        }

        for( int index = 0; index < sher_mem.variants.size(); ++index )
        {
            bool null(false), splice(false), missense(false), initiator(false),
                silent(false), intronic(false), inframe(false), noncoding(false); 
            for(auto& var_type: sher_mem.variants[index].type)
            {
                if(null_types.contains(var_type))
                    null = true;
                if(splice_types.contains(var_type))
                    splice = true;
                if(missense_types.contains(var_type))
                    missense = true;
                if(initiator_types.contains(var_type))
                    initiator = true;
                if(silent_types.contains(var_type))
                    silent = true;
                if(intronic_types.contains(var_type))
                    intronic = true;
                if(inframe_types.contains(var_type))
                    inframe = true;
                if(noncoding_types.contains(var_type))
                    noncoding = true;
            }
            // I think variant might go multiple variant type subtrees
            // e.g. "splice_region_variant;splice_polypyrimidine_tract_variant;intron_variant;non_coding_transcript_variant"
            // this kinda variant should go splice to check if it disrupt donor or acceptor,
            // but it can also go to intronic to gain some benign score.
            // FIXME: Am I wrong?
            if(null)
                go_null( sher_mem, index, db.db_gtf, para, db, sher_conseq );
            if(splice)
                go_splice( sher_mem, index, db.db_gtf, para, db.db_fasta );
            if(initiator)
                go_initiator( sher_mem, index, para, db.db_uniprot );
            if(missense)
                go_missense( sher_mem, index, para, db.db_uniprot, db.db_gtf, db.db_clinvar);
            if(silent)
                go_silent( sher_mem, index, db.db_gtf );
            if(intronic)
                go_intronic( sher_mem, index, db.db_gtf );
            if(inframe)
                go_inframe( sher_mem, index, para, sher_conseq);
            if(noncoding)
                // FIXME: in sherloc-table, EV0168:
                //      "... Do not use for structural or functional RNAs 
                //      that do not code for any known protein (i.e. TERC RNA)."
                // So maybe not all types of variant not in the above list are going to non_coding
                go_non_coding( sher_mem, index, para );
        }
    }

    void run( Patient::Patient& ze
                  , DB::DBSet& db  
                  , SherlocConsequence& sher_conseq 
                )
    {
        decltype(auto) para = SherlocParameter::get_paras();
        for ( auto& sher_mem: ze.sher_mems ) {
            run( sher_mem, para, db, sher_conseq );
        }
    }
    
};

}
