/**
 * @file option_parser.hpp
 * @brief The program option parser definition.
 */
#pragma once
#include <boost/program_options.hpp>
#include <stdexcept>
namespace nucleona::app::cli{

namespace po = boost::program_options;
/**
 * @brief The option parser base type.
 * @details Inheritance this OptionParser and customize option parser for client application.
 */
class OptionParser
{
protected:
    po::variables_map vm;
public:
    /**
     * @brief Get the command line parameter by parameter name.
     * 
     * @param parameter_name Name of the parameter.
     * @param variable The variable will be assigned.
     * @tparam  StringType Type of string, will be auto reduced.
     * @tparam ParType  Parameter type of variable, will be auto reduced. 
     */
    template< class StringType, class ParType >
    void get_parameter(StringType&& parameter_name, ParType& variable)
    {
        if(vm.count(parameter_name))
            variable = vm[parameter_name].template as<ParType>();
    }
    /**
     * @brief Get the command line parameter by parameter name.
     * @details Boolean flag version.
     * 
     * @param parameter_name Name of the parameter.
     * @param variable The variable will be assigned.
     * @tparam  StringType  Type of string, will be auto reduced.
     */
    template< class StringType >
    void get_parameter(StringType&& parameter_name, bool& variable)
    {
        variable = vm.count(parameter_name);
    }
    /**
     * @brief Set the parameter which must input.
     * @details Some parameter is necessary for program. This method is use to verify the parameters which is necessary.
     * 
     * @param var_name The variable which is necessary. 
     * @param x The Variable which stored the parameter.
     * @param x_def The variable initial status.
     * @tparam  T1 Variable type, will be reduced.
     * @tparam String  Variable name type ( usually string ), will be reduced.
     */
    template< class T1, class T2, class String >
    void necessary_verification(String&& var_name, T1&& x, T2&& x_def )
    {
        if( x == x_def )
        {
            throw std::runtime_error("error: option " + var_name + " must be given a parameter.");
        }
    }
};
}