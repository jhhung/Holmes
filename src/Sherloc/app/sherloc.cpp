#include <iostream> 
#include <Sherloc/app/sherloc/main.hpp>

int main( int argc, const char* argv[] )
{
    Sherloc::app::sherloc::GetParameters parameters( argc, argv );    
    Sherloc::app::sherloc::Run( parameters );
    return 0;
}
