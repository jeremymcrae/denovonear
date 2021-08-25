#ifndef GTF_H_
#define GTF_H_

#include <cstdint>
#include <fstream>
#include <string>

#include "gzstream/gzstream.h"

namespace gencode {

// store required fields from a GTF line
struct GTFLine {
    std::string chrom;
    std::string feature;
    int start;
    int end;
    std::string strand;
    std::string symbol;
    std::string tx_id;
    std::string transcript_type;
    bool is_principal;
};

GTFLine parse_gtfline(std::string &line);

class GTF
{
    // we need to allow for gzipped GTFs or not, so prepare handles for both.
    // This is wasteful, but only one will open and unsure of an easier way.
    std::ifstream handle;
    igzstream gzhandle;  
    std::string line;
    bool gzipped;
public:
    GTF(std::string path);
    GTFLine next();
};

} // namespace

#endif // GTF_H_