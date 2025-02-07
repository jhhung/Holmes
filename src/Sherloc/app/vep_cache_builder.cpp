#include <iostream> 
#include <Sherloc/app/vep_cache_builder/main.hpp>

int main( int argc, const char* argv[] )
{
    Sherloc::app::vep_cache_builder::GetParameters parameters( argc, argv );    
    Sherloc::app::vep_cache_builder::Run( parameters );
    return 0;
}
