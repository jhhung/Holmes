#include <iostream> 
#include <Sherloc/app/disease_json_builder/main.hpp>

int main( int argc, const char* argv[] )
{
    Sherloc::app::disease_json_builder::GetParameters parameters( argc, argv );    
    Sherloc::app::disease_json_builder::Run( parameters );
    return 0;
}
