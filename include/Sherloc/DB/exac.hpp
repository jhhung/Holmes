#pragma once

#include <vector>
#include <string>

namespace Sherloc::DB {

class Exac
{
  public:
    std::vector<Exac> other;

    size_t chr; 

    size_t pos;

    std::string ref;
    
    std::string alt;

    char status;

    char hom;

    Exac():chr(0), pos(0), ref(), alt(), status('0'), hom('0'){}

    Exac( const std::string& ref0, const std::string& alt0 )
        :chr(0), pos(0), ref(ref0), alt(alt0), status('0'), hom('0'){}

    Exac& operator = ( const Exac& a )
    {
        chr = a.chr;
        pos = a.pos;
        ref = a.ref;
        alt = a.alt;
        status = a.status;
        hom = a.hom;
        return *this;
    }

    void clear()
    {
        Exac a;
        *this = a;
    }
};

}
