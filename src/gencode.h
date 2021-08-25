#ifndef GENCODE_H_
#define GENCODE_H_

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include <set>

#include "gtf.h"
#include "tx.h"

namespace gencode {

// collects info for a transcript, so we can later construct a Tx object
struct TxInfo {
    std::string name = "";
    std::string chrom;
    int start;
    int end;
    std::string strand;
    std::vector<std::vector<int> > exons;
    std::vector<std::vector<int> > cds;
    int offset = 0;
};

// stores HGNC symbol with the transcript, so we can collect transcripts by gene
struct NamedTx {
    std::string symbol;
    Tx tx;
    bool is_principal;
};

// open the gencode annotations GTF, and return Tx objects for each transcript
std::vector<NamedTx> open_gencode(std::string path, bool coding=true);

} // namespace

#endif // GENCODE_H_

