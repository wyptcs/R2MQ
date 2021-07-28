#ifndef BPINDEXING_HPP
#define BPINDEXING_HPP
#include "sdsl_modified/bp_support_sada.hpp"
#include "sdsl_modified/rank_support_v5.hpp"
#include <sdsl/rrr_vector.hpp>
#include <sdsl/rmq_succinct_sada.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdint.h>
using namespace std;
using namespace sdsl;
class BPIndexing {
    typedef int arraysize;
    static const int BLOCK_SIZE = 512;
    arraysize n;

public:
    bit_vector bp;
    bp_support_sada<BLOCK_SIZE,1,rank_support_v5<>,select_support_mcl<1,1>> bps;
    BPIndexing(vector<arraysize>& array);
    double getSpaceRMQ();
    void makeBP(vector<arraysize>& array);
    arraysize rmq(arraysize i,arraysize j);
    void initializeCache();
};
#endif