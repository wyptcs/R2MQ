#ifndef DFUDSINDEXING_HPP
#define DFUDSINDEXING_HPP
#include "sdsl_modified/bp_support_sada.hpp"
#include "sdsl_modified/rank_support_v5.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdint.h>
using namespace std;
using namespace sdsl;
class DFUDSIndexing {
    typedef uint arraysize;
    static const int BLOCK_SIZE = 512;
    arraysize n;
private:
    void makeDFUDS(vector<arraysize>& array);
    arraysize rmqe(arraysize i,arraysize j);
public:
    bit_vector dfuds;
    bp_support_sada<BLOCK_SIZE,1,rank_support_v5<>,select_support_mcl<0,1>> bps;
    DFUDSIndexing(vector<arraysize>& array);
    arraysize rmq(arraysize i,arraysize j);
    double getSpaceRMQ();
    void initializeCache();
};
#endif