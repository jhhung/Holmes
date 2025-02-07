#pragma once
#include <iostream>
#include <fstream>
#include <Sherloc/option_parser.hpp>
#include <nlohmann/json.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace Sherloc {
namespace app {
namespace disease_json_builder {

struct Parameters
{
    int id;
    std::string file;
    bool adar;
    bool onset;
    bool severe;
    double yield;
    std::string output;
};

class GetParameters :
    public Parameters,
    public nucleona::app::cli::OptionParser
{
public:
    GetParameters( int argc, char const * argv[] )
    {
        adar = false;
        onset = false;
        severe = false;
        namespace po = boost::program_options;
        po::options_description desc( "Allowed options") ;
        desc.add_options()
            ( "help,h"  , "show help message" )
            ( "id,i", po::value<int >()->required(), "disease id" )
            ( "file,f", po::value< std::string >()->required(), "disease list archive file" )
            ( "adar,a"   , po::bool_switch( &adar ), "is it AD  ?" )
            ( "onset,n"   , po::bool_switch( &onset ), "is it EARLY ?" )
            ( "severe,s"   , po::bool_switch( &severe), "is it severe ?" )
            ( "yield,y"   , po::value< double >()->default_value(0.0), "yield test ( ? % )" )
            ( "output,o"   , po::value< std::string >()->default_value("/tmp/disease.json"), "output data" )
        ;

        po::store( po::parse_command_line( argc, argv, desc ), vm );

        if( argc == 1 or vm.count( "help" ))
        {
            std::cout << desc << std::endl;
            std::exit(1);
        }

        po::notify( vm );

        get_parameter( "id", id );
        get_parameter( "file", file );
        get_parameter( "yield", yield );
        get_parameter( "output" , output );
    }
};

template< class GET_PARAMETERS >
void Run( GET_PARAMETERS&& args )
{
    nlohmann::json disease;

    std::ofstream output( args.output );

    std::cout << "Create json object ...\n";

    disease[ "id" ]  = args.id;
    disease[ "file" ] = args.file;
    disease[ "adar" ] = args.adar;
    disease[ "onset" ] = args.onset;
    disease[ "severe" ] = args.severe;
    disease[ "yield" ] = args.yield;

    std::cout << "Output json object to " << args.output << '\n';

    output << disease.dump(4);
    output.close();
}

}}}
