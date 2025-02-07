#include <Sherloc/app/archive_compressor/main.hpp>

int main( int argc, const char* argv[] )
{
    Sherloc::app::archive_compressor::GetParameters parameters( argc, argv );    
    Sherloc::app::archive_compressor::Run( std::move(parameters) );
    return 0;
}
