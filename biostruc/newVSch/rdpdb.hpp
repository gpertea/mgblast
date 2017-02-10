#ifndef _RDPDB_HPP
#define _RDPDB_HPP

#include <cgi/cgictx.hpp>

using namespace ncbi;
using namespace std;


namespace RdPdb {
	
    const int ATMNUM_WIDTH = 5;
    const int ATMNAME_WIDTH = 4;
    const int RESNAME_WIDTH = 4;
    const int RESNUM_WIDTH = 4;
    const int COORD_WIDTH = 8;
    const int OCC_WIDTH = 6;
    const int BVAL_WIDTH = 6;
    const int SEQNUM_WIDTH = 4;
    const int CYSNAME_WIDTH = 3;
    const int CHAIN_WIDTH = 1;
    const int BUFFSIZE = 1024;
    const int FIRST_WIDTH = 6;
};


bool ReadPDB(CCgiContext& ctx);

#endif
