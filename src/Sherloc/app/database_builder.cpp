#include <iostream> 
#include <Sherloc/app/database_builder/main.hpp>

int main( int argc, const char* argv[] )
{
    Sherloc::app::database_builder::GetParameters parameters( argc, argv );    
    Sherloc::app::database_builder::Run( parameters );
    return 0;
}
