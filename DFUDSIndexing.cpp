#include "DFUDSIndexing.hpp"
#include "sdsl_modified/bp_support_sada.hpp"
#include "sdsl_modified/rank_support_v5.hpp"
#include <iostream>
#include <vector>
#include <stack>
using namespace std;
using namespace sdsl;
DFUDSIndexing::DFUDSIndexing(vector<DFUDSIndexing::arraysize> &array) {
    n = array.size();
    makeDFUDS(array);
    cout<<"DFUDS BUILD FINISHED"<<endl;
    n = array.size();
    bp_support_sada<DFUDSIndexing::BLOCK_SIZE,1,rank_support_v5<>,select_support_mcl<0,1>> t(&dfuds);
    bps = t;
    cout<<"RMM TREE BUILD FINISHED"<<endl;
}
void DFUDSIndexing::makeDFUDS(vector<DFUDSIndexing::arraysize> & array) {
    dfuds.resize(array.size()*2+2);
    DFUDSIndexing::arraysize pos = array.size()*2+1;
    stack<DFUDSIndexing::arraysize> Q;
    for(DFUDSIndexing::arraysize i = array.size()-1;;i--){
        dfuds[pos] = 0;
        pos--;
        while((not Q.empty()) and Q.top() <= array[i]){
            dfuds[pos] = 1;
            pos--;
            Q.pop();
        }
        Q.push(array[i]);
        if(i==0)
            break;
    }
    dfuds[pos] = 0;
    pos--;
    while(not Q.empty()){
        dfuds[pos] = 1;
        pos--;
        Q.pop();
    }
    dfuds[pos] = 1;
}
double DFUDSIndexing::getSpaceRMQ() {
    int64_t dfuds_size = size_in_bytes(dfuds)*8;
    int64_t bps_size = size_in_bytes(bps)*8;
    return sizeof(DFUDSIndexing::arraysize)*8+double(dfuds_size + bps_size)/double(n);
}
DFUDSIndexing::arraysize DFUDSIndexing::rmqe(DFUDSIndexing::arraysize i, DFUDSIndexing::arraysize j) {//i<=x<=j
    if(i==j)
        return bps.select(i+1);
    DFUDSIndexing::arraysize l = bps.select(i+2);
    DFUDSIndexing::arraysize r = bps.select(j+1);
    return bps.rmq_left(l,r);
}
DFUDSIndexing::arraysize DFUDSIndexing::rmq(DFUDSIndexing::arraysize i, DFUDSIndexing::arraysize j) {//return plus 1 for start at 0
    if(i==j)
        return i;
    DFUDSIndexing::arraysize w = rmqe(i,j);
    if(bps.rank0(bps.find_open(w)) == i+1){
        return i;
    }
    return bps.rank0(w)-1;
}
void DFUDSIndexing::initializeCache() {
    bps.initializeCache();
}