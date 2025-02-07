#include <iostream> 
#include <Sherloc/app/vep_cache_query/main.hpp>

int main( int argc, const char* argv[] )
{
    Sherloc::app::vep_cache_query::GetParameters parameters( argc, argv );    
    Sherloc::app::vep_cache_query::Run( parameters );
    return 0;
}
