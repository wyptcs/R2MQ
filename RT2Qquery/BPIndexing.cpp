#include "BPIndexing.hpp"
#include "sdsl_modified/bp_support_sada.hpp"
#include "sdsl_modified/rank_support_v5.hpp"
#include <sdsl/rrr_vector.hpp>
#include <sdsl/rmq_succinct_sada.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <stack>
#include <cassert>
#include <math.h>
#include <unordered_map>
using namespace std;
using namespace sdsl;
BPIndexing::BPIndexing(vector<BPIndexing::arraysize> &array) {
    n = array.size();
    makeBP(array);
    cout<<"BP BUILD FINISHED"<<endl;
    bp_support_sada<BPIndexing::BLOCK_SIZE,1,rank_support_v5<>,select_support_mcl<1,1>> t(&bp);
    bps = t;
    cout<<"RMM TREE BUILD FINISHED"<<endl;
}
void BPIndexing::makeBP(vector<BPIndexing::arraysize> & array) {
    bp.resize(array.size()*2+2);
    BPIndexing::arraysize pos = 0;
    stack<BPIndexing::arraysize> Q;
    bp[pos] = 1;
    pos++;
    bp[pos] = 1;
    pos++;
    Q.push(array[0]);
    for(BPIndexing::arraysize i = 1;i < array.size();i++){
        while((not Q.empty()) and Q.top() < array[i]){
            bp[pos]=0;
            pos++;
            Q.pop();
        }
        Q.push(array[i]);
        bp[pos]=1;
        pos++;
    }
    while(not Q.empty()){
        bp[pos] = 0;
        pos++;
        Q.pop();
    }
    bp[pos]=0;
    pos++;
#ifdef DEBUG
    assert(pos==array.size()*2+2);
#endif
}

BPIndexing::arraysize BPIndexing::rmq(BPIndexing::arraysize i, BPIndexing::arraysize j) {
    return bps.rank(bps.rmq(bps.select(i + 2) - 1, bps.select(j + 2)))-1;
}
double BPIndexing::getSpaceRMQ() {
    int64_t bp_size = size_in_bytes(bp)*8;
    int64_t bps_size = size_in_bytes(bps)*8;
    return sizeof(BPIndexing::arraysize)*8 + double(bp_size + bps_size)/double(n);
}
void BPIndexing::initializeCache() {
    bps.initializeCache();
}
