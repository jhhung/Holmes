#pragma once

#include <vector>
#include <string>
#include <nlohmann/json.hpp>
#include <Sherloc/option_parser.hpp>
#include <Sherloc/DB/db.hpp>


namespace Sherloc {

class Position
{
  public:

    std::string chr;

    int start;

    int end;

    std::string gene;

    std::string id;

    Position(){};

    Position ( const std::string& chr0, const int& start0, const int& end0, const std::string& gene0, const std::string& id0 ):
        chr( chr0 ), start( start0 ), end( end0 ), gene( gene0 ), id( id0 ){}
};

class Gene_pos
{
  public:
    std::string chr;

    size_t start;
    
    size_t end;

    bool adar;

    bool onset;

    bool severe;

    bool yield;
};

class Disease 
{
  public:

    bool empty;

    int id;

    std::vector< Position > pos;

    std::vector< Gene_pos > gene_pos;

    bool adar;

    bool onset;

    bool severe;

    double yield;

    Disease():empty(true), adar(false), onset(false), severe(false), yield(0.0){}

    Disease( const nlohmann::json& js, const std::filesystem::path& disease_db_file)
        :empty(false), adar(false), onset(false), severe(false), yield(0.0)
    {
        id = js["id"];
        adar = js["adar"];
        onset = js["onset"];
        severe = js["severe"];
        yield = js["yield"];

        if( id == 0 )
        {
            empty = true;
            return;
        }

        std::map<int, std::vector< std::map< std::string, std::string > > > omim;
        DB::load_archive_from(omim, disease_db_file);
        
        auto it = omim.find( id );
        if( it == omim.end() )  *this = Disease();
        else
        {
            auto vec = it->second;
            for( auto& i: vec )
            {
                /* if vcf file chrosome column have prefix "chr" */
                std::string chr, chr_no;
                chr = i.find("chr")->second;
                if (chr.substr(0, 3) == "chr") {
                  chr_no = std::move(chr.substr(3));
                } else {
                  chr_no = chr;
                }
                // Unknwon chromosome or mitochondrial
                if (chr_no[0] == 'U' or chr_no[0] == 'M') {
                  continue;
                }
                pos.emplace_back( Position( chr_no, std::stoi( i.find("start")->second ), std::stoi( i.find("end")->second ), i.find("gene")->second, i.find("id")->second ) );
            }
        }

        for( auto& i: js["gene"] )
        {
            Gene_pos a;
            a.chr = i["chr"];
            a.start = i["start"];
            a.end = i["end"];
            a.adar = i["adar"];
            a.onset = i["onset"];
            a.severe = i["severe"];
            a.yield = i["yield"];
            gene_pos.emplace_back( a );
        }
    }
};

}
