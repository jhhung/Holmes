#pragma once

#include <iostream>
#include <fstream>
#include <filesystem>
#include <Sherloc/option_parser.hpp>
#include <Sherloc/DB/dbset.hpp>
#include <spdlog/spdlog.h>

namespace Sherloc::app::database_builder {

struct Parameters
{
    std::string clinvar_file;
    std::string k_file;
    std::string coverage_file;
    std::pair<std::string, std::string> gene_info_file;
    std::string dvd_file;
    std::string uniprot_file;
    std::string fasta_file;
    std::string ensembl_gtf_file;
    std::string refseq_gtf_file;
    std::string gnom_file;
    std::string output;
    bool separate;
    int thread_num;
};

class GetParameters :
    public Parameters,
    public nucleona::app::cli::OptionParser
{
public:
    GetParameters( int argc, char const * argv[] )
    {
        separate = false;
        namespace po = boost::program_options;
        po::options_description desc( "Allowed options") ;
        desc.add_options()
            ( "help,h",     "show help message" )

            ( "k,k",        po::value<std::string>(&k_file)->default_value(""), "1000 Genomes Project file" )

            ( "clinvar,c",  po::value<std::string>(&clinvar_file)->default_value(""), "VEP annotated Clinvar file" )

            ( "coverage,v", po::value<std::string>(&coverage_file)->default_value(""), "Coverage file" )

            ( "omim",       po::value<std::string>(&gene_info_file.first)->default_value(""), "Gene info: OMIM genemap file" )
            ( "ncbigene",   po::value<std::string>(&gene_info_file.second)->default_value(""), "Gene info: NCBI gene file" )

            ( "dvd,d",      po::value<std::string>(&dvd_file)->default_value(""), "DVD database file" )

            ( "uniprot,u",  po::value<std::string>(&uniprot_file)->default_value(""), "UniProt database file" )

            ( "fasta,f",    po::value<std::string>(&fasta_file)->default_value(""), "Reference Fasta file" )

            ( "ensembl_gtf",po::value<std::string>(&ensembl_gtf_file)->default_value(""), "Ensembl GTF file" )
            ( "refseq_gtf", po::value<std::string>(&refseq_gtf_file)->default_value(""), "Refseq GTF file" )

            ( "gnom_file,g", po::value<std::string>(&gnom_file)->default_value(""),
                "gnomAD url list file (each chr url/path per line)" )
            
            ( "output,o",   po::value<std::string>(&output)->default_value("/tmp/"), "Output directory" )

            ( "thread,t",   po::value<int>(&thread_num)->default_value(4), "Thread num for parallel building GnomAD" )
        ;

        po::store( po::parse_command_line( argc, argv, desc ), vm );

        if( argc == 1 or vm.count( "help" ))
        {
            std::cout << desc << std::endl;
            std::exit(1);
        }

        po::notify( vm );
    }
};


template< class GET_PARAMETERS >
void Run( GET_PARAMETERS&& args )
{
    // 1. Read database in
    // 2. Create database object
    // 3. Save database objects into one archive file
    using namespace std::filesystem;
    Sherloc::DB::DBSet dbset;
    std::string database_to_build = "[";
    auto output_dir = path(args.output);
    if(!exists(output_dir)){
        SPDLOG_ERROR("Output dir: {} not exists!", output_dir.string());
        exit(1);
    }

    if(!args.k_file.empty()){
        database_to_build += "K Genome,";
        SPDLOG_INFO("Building K Genome...");
        dbset.db_1kg.from(args.k_file);
        dbset.db_1kg.save(output_dir / DB::DBSet::get_default("1kg"));
    }
    
    if(!args.clinvar_file.empty()){
        database_to_build += "ClinVar,";
        SPDLOG_INFO("Building ClinVar...");
        dbset.db_clinvar.from(args.clinvar_file);
        dbset.db_clinvar.save(output_dir / DB::DBSet::get_default("clinvar"));
    }

    if(!args.gene_info_file.first.empty() and !args.gene_info_file.first.empty()){
        database_to_build += "Gene Info,";
        SPDLOG_INFO("Building Gene Info...");
        nlohmann::json gene_info_config = {
            {"omim", args.gene_info_file.first},
            {"ncbi", args.gene_info_file.second}
        };
        auto tmp_config_file = std::filesystem::temp_directory_path() / "holmes_gene_info.json";
        {
            auto os = std::ofstream{tmp_config_file};
            os << gene_info_config;
        }
        dbset.db_gene_info.from(tmp_config_file);
        dbset.db_gene_info.save(output_dir / DB::DBSet::get_default("gene_info"));
    }

    if(!args.dvd_file.empty()){
        database_to_build += "DVD,";
        SPDLOG_INFO("Building DVD...");
        dbset.db_dvd.from(args.dvd_file);
        dbset.db_dvd.save(output_dir / DB::DBSet::get_default("dvd"));
    }

    if(!args.uniprot_file.empty()){
        database_to_build += "UniProt,";
        SPDLOG_INFO("Building UniProt...");
        dbset.db_uniprot.from(args.uniprot_file);
        dbset.db_uniprot.save(output_dir / DB::DBSet::get_default("uniprot"));
    }

    if(!args.fasta_file.empty()){
        SPDLOG_WARN("Fasta builder is deprecated. straight up load indexed fasta.gz instead.");
        // database_to_build += "Fasta,";
        // SPDLOG_INFO("Building Fasta...");
        // dbset.db_fasta.from(args.fasta_file);
        // dbset.db_fasta.save(output_dir / DB::DBSet::get_default("fasta"));
    }

    if(!args.ensembl_gtf_file.empty() and !args.refseq_gtf_file.empty()){
        database_to_build += "GTF,";
        SPDLOG_INFO("Building GTF...");
        nlohmann::json gtf_config = {
            {"ensembl", args.ensembl_gtf_file},
            {"refseq", args.refseq_gtf_file}
        };
        auto tmp_config_file = std::filesystem::temp_directory_path() / "holmes_gtf.json";
        {
            auto os = std::ofstream{tmp_config_file};
            os << gtf_config;
        }
        dbset.db_gtf.from(tmp_config_file);
        dbset.db_gtf.save(output_dir / DB::DBSet::get_default("gtf"));
    }

    if(!args.coverage_file.empty()){
        database_to_build += "Coverage,";
        SPDLOG_INFO("Building Coverage...");
        dbset.db_coverage.from(args.coverage_file);
        dbset.db_coverage.save(output_dir / DB::DBSet::get_default("coverage"));
    }

    if(!args.gnom_file.empty()){
        database_to_build += "GnomAD,";
        SPDLOG_INFO("Building GnomAD...");
        dbset.db_gnom.gnom_dir = output_dir / DB::DBSet::get_default("gnom");
        dbset.db_gnom.thread_num = std::max(1, args.thread_num);
        dbset.db_gnom.from(args.gnom_file);
    }
    
    if(database_to_build.back() == ',')
        database_to_build.pop_back();
    database_to_build += "]";
    SPDLOG_INFO("Built Databases: {}", database_to_build);
    SPDLOG_INFO("Archive files go to: {}", output_dir.string());
}

}
