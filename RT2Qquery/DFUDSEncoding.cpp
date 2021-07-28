#include "DFUDSEncoding.hpp"
#include "sdsl_modified/bp_support_sada.hpp"
#include "sdsl_modified/rank_support_v5.hpp"
#include <sdsl/rrr_vector.hpp>
#include <sdsl/rmq_succinct_sada.hpp>

#include <iostream>
#include <vector>
#include <algorithm>
#include <stack>
#include <cassert>
#include <unordered_map>
#include <math.h>
#include <cassert>
using namespace std;
using namespace sdsl;
int INITSTART = 5;
int INITINTERVAL = 10;
int INITQUERYSIZE = 10000;
MODE getMode(int value){
    if(value <= UINT8_MAX)
        return CHAR;
    else if(value <= UINT16_MAX)
        return SHORT;
    else
        return INT;
}
DFUDSEncoding::DFUDSEncoding(vector<DFUDSEncoding::arraysize> &array,double usingSpace,bool reverse) {
    this->reverse = reverse;
    this->isSavingDepth = false;
    makeRMQStructure(array);
    makeTop2Structure(array);
    makeDepthPos(array.size(), INITSTART, INITINTERVAL);
    //cout << "DEPTH BUILD FINISHED" << endl;
    makeLDepthPos(array.size(), INITSTART, INITINTERVAL);
    //cout << "LDEPTH BUILD FINISHED" << endl;
    allDepthValues.clear();
    maxDepth = *max_element(depthValues.begin(),depthValues.end());
    int maxLDepth = *max_element(lDepthValues.begin(),lDepthValues.end());
    this->interval = 0;
    cout<<"maxdepth :"<<maxDepth<<endl;
    if(maxDepth > log2(array.size()) and usingSpace > 0) {
        this->isSavingDepth = true;
        MODE depthMode = getMode(maxDepth);
        MODE lDepthMode = getMode(maxLDepth);
        GetInterval getInterval(array, int(depthMode), int(lDepthMode),INITQUERYSIZE);
        auto startPosAndInterval = getInterval.getStartPosAndInterval(usingSpace);
        depthValues.clear();
        lDepthValues.clear();

        /*makeDepthPos(array.size(), startPosAndInterval.first, startPosAndInterval.second);
        makeLDepthPos(array.size(), startPosAndInterval.first, startPosAndInterval.second);*/
        //vector<int> tempSaving = {1,2,3,4,5,6,7,8,9,10,20};
        vector<int> saving;
        saving = getInterval.getSavingCandidateOne(usingSpace);
        makeDepthPosNew(array.size(),saving);
        makeLDepthPosNew(array.size(),saving);

        this->interval = startPosAndInterval.second;
        allDepthValues.clear();
        makeDepthValuesContainer(depthMode);
        makeLDepthValuesContainer(lDepthMode);
        depthValues.clear();
        lDepthValues.clear();
        //cout << "DFUDS TOP2 STRUCTURE WITH SAVING DEPTH FINISHED" << endl;
        //cout<<"starting point : "<<startPosAndInterval.first<<endl;
        //cout<<"interval : "<<startPosAndInterval.second<<endl;
    }else{
        bit_vector r(1,0);
        depthPos = rrr_vector<>(r);
    }
}
DFUDSEncoding::DFUDSEncoding(vector<arraysize> &array) {
    this->isSavingDepth = false;
    makeRMQStructure(array);
    makeTop2Structure(array);
}
DFUDSEncoding::DFUDSEncoding() {

}

void DFUDSEncoding::makeDepthValuesContainer(MODE mode) {
    depthValuesContainer = Depth(mode);
    for(arraysize i = 0;i<depthValues.size();i++)
        depthValuesContainer.push(depthValues[i]);
}
void DFUDSEncoding::makeLDepthValuesContainer(MODE mode) {
    lDepthValuesContainer = Depth(mode);
    for(arraysize i = 0;i<lDepthValues.size();i++)
        lDepthValuesContainer.push(lDepthValues[i]);
}
void DFUDSEncoding::makeRMQStructure(vector<DFUDSEncoding::arraysize> &array) {
    makeDFUDS(array);
    //cout<<"DFUDS BUILD FINISHED"<<endl;
    n = array.size();
    bp_support_sada<DFUDSEncoding::BLOCK_SIZE,1,rank_support_v5<>,select_support_mcl<0,1>> t(&dfuds);
    bps = t;
    //cout<<"RMM TREE BUILD FINISHED"<<endl;
}
void DFUDSEncoding::makeTop2Structure(vector<DFUDSEncoding::arraysize> &array){
    lLeavesRank = rank_support_v5<100,3>(&dfuds);
    //cout<<"LLEAVESRANK BUILD FINISHED"<<endl;
    makeSpineUsingTwoDVector(array);
    spineRank = rank_support_v5<>(&spine);
    //cout<<"SPINE BUILD FINISHED"<<endl;
    spine_ss_0 = select_support_mcl<0,1>(&spine);
    spine_ss_1 = select_support_mcl<1,1>(&spine);
    //cout<<"SPINE SELECT BUILD FINISHED"<<endl;
}
void DFUDSEncoding::makeDepthPosNew(DFUDSEncoding::arraysize n,vector<int> savePoint) {
    bit_vector depthPosRaw(bps.size(),0);
    isCompleteDepth.resize(bps.size());
    allDepthValues.push_back(0);
    for(DFUDSEncoding::arraysize i = 0;i<n;i++){
        DFUDSEncoding::arraysize treeIndex = arrayIndexToTreeIndex(i+1);
        int depth = depth_init(treeIndex);
        if(find(savePoint.begin(),savePoint.end(),depth) != savePoint.end()){
            depthValues.push_back(depth);
            depthPosRaw[treeIndex] = 1;
        }
    }
    isCompleteDepth.resize(0);
    depthPos = rrr_vector<>(depthPosRaw);
    depthPosRank = rrr_vector<>::rank_1_type(&depthPos);
}
void DFUDSEncoding::makeLDepthPosNew(DFUDSEncoding::arraysize n,vector<int> savePoint){
    for(DFUDSEncoding::arraysize i = 0;i<n;i++){
        DFUDSEncoding::arraysize treeIndex = arrayIndexToTreeIndex(i+1);
        int node_depth = allDepthValues[i+1];
        if(find(savePoint.begin(),savePoint.end(),node_depth) != savePoint.end()){
            lDepthValues.push_back(getLdepth_init(treeIndex));
        }
    }
}
void DFUDSEncoding::makeDepthPos(DFUDSEncoding::arraysize n,int startDepth,int interval) {
    bit_vector depthPosRaw(bps.size(),0);
    isCompleteDepth.resize(bps.size());
    allDepthValues.push_back(0);
    //cout<<startDepth<<endl;
    for(DFUDSEncoding::arraysize i = 0;i<n;i++){
        DFUDSEncoding::arraysize treeIndex = arrayIndexToTreeIndex(i+1);
        int depth = depth_init(treeIndex);
        if(depth >= startDepth and depth%interval == startDepth%interval){
        //if((startDepth !=INITSTART and depth <= startDepth) or (startDepth == INITSTART and depth >= startDepth and depth%interval == startDepth%interval)){
            depthValues.push_back(depth);
            depthPosRaw[treeIndex] = 1;
        }
    }
    isCompleteDepth.resize(0);
    depthPos = rrr_vector<>(depthPosRaw);
    depthPosRank = rrr_vector<>::rank_1_type(&depthPos);
}
void DFUDSEncoding::makeLDepthPos(DFUDSEncoding::arraysize n, int startDepth, int interval) {
    for(DFUDSEncoding::arraysize i = 0;i<n;i++){
        DFUDSEncoding::arraysize treeIndex = arrayIndexToTreeIndex(i+1);
        int node_depth = allDepthValues[i+1];
        if(node_depth >= startDepth and node_depth%interval == startDepth%interval){
        //if((startDepth !=INITSTART and node_depth <= startDepth) or (startDepth == INITSTART and node_depth >= startDepth and node_depth%interval == startDepth%interval)){
            lDepthValues.push_back(getLdepth_init(treeIndex));
        }
    }
}
void DFUDSEncoding::makeDFUDS(vector<DFUDSEncoding::arraysize> & array) {
    dfuds.resize(array.size()*2+2);
    DFUDSEncoding::arraysize pos = array.size()*2+1;
    stack<DFUDSEncoding::arraysize> Q;
    for(DFUDSEncoding::arraysize i = array.size()-1;;i--){
        dfuds[pos] = 0;
        pos--;
        if(reverse){
            while((not Q.empty()) and Q.top() < array[i]){
                dfuds[pos] = 1;
                pos--;
                Q.pop();
            }
        }else {
            while ((not Q.empty()) and Q.top() <= array[i]) {
                dfuds[pos] = 1;
                pos--;
                Q.pop();
            }
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
void DFUDSEncoding::makeSpineUsingTwoDVector(vector<DFUDSEncoding::arraysize> &array) {
    vector<vector<DFUDSEncoding::arraysize>> leftSpineTable;
    vector<vector<DFUDSEncoding::arraysize>> rightSpineTable;
    for(DFUDSEncoding::arraysize i = 0;i<array.size();i++){
        vector<DFUDSEncoding::arraysize> p;
        leftSpineTable.push_back(p);
        rightSpineTable.push_back(p);
    }
    stack<DFUDSEncoding::arraysize> Q;
    for(DFUDSEncoding::arraysize i = 0;i < array.size();i++) {
        if(reverse){
            while (not(Q.empty()) and array[Q.top()] <= array[i]) {
                leftSpineTable[i].push_back(Q.top());
                Q.pop();
            }
            if(not(Q.empty()) and array[Q.top()] > array[i]) {
                rightSpineTable[Q.top()].push_back(i);
            }
        }else {
            while (not(Q.empty()) and array[Q.top()] < array[i]) {
                leftSpineTable[i].push_back(Q.top());
                Q.pop();
            }
            if (not(Q.empty()) and array[Q.top()] >= array[i]) {
                rightSpineTable[Q.top()].push_back(i);
            }
        }
        Q.push(i);
    }
    DFUDSEncoding::arraysize spine_size = 0;
    for(DFUDSEncoding::arraysize i = 0;i < array.size();i++){
        int64_t leftsize = leftSpineTable[i].size();
        int64_t rightsize = rightSpineTable[i].size();
        spine_size += max(leftsize + rightsize - 1, (int64_t)0);
    }
    spine.resize(spine_size+2);
    DFUDSEncoding::arraysize pos = 0;
    for(DFUDSEncoding::arraysize i = 0;i < array.size();i++){
        int64_t leftsize = leftSpineTable[i].size();
        int64_t rightsize = rightSpineTable[i].size();
        int64_t left_index = leftSpineTable[i].size() - 1;
        int64_t right_index = rightSpineTable[i].size() - 1;
        if(reverse)
        {
            for(int64_t j = 0;j<leftsize + rightsize - 1;j++){
                if(left_index == -1){
                    spine[pos] = 1;
                    pos++;
                    right_index--;
                }else if(right_index == -1) {
                    spine[pos] = 0;
                    pos++;
                    left_index--;
                }else if(array[leftSpineTable[i][left_index]] <= array[rightSpineTable[i][right_index]]){
                    spine[pos] = 1;
                    pos++;
                    right_index--;
                }else if(array[leftSpineTable[i][left_index]] > array[rightSpineTable[i][right_index]]){
                    spine[pos] = 0;
                    pos++;
                    left_index--;
                }else{
                    cout<<"ERROR ON BUILD SPINE"<<endl;
                }
            }
        }else {
            for (int64_t j = 0; j < leftsize + rightsize - 1; j++) {
                if (left_index == -1) {
                    spine[pos] = 1;
                    pos++;
                    right_index--;
                } else if (right_index == -1) {
                    spine[pos] = 0;
                    pos++;
                    left_index--;
                } else if (array[leftSpineTable[i][left_index]] < array[rightSpineTable[i][right_index]]) {
                    spine[pos] = 1;
                    pos++;
                    right_index--;
                } else if (array[leftSpineTable[i][left_index]] >= array[rightSpineTable[i][right_index]]) {
                    spine[pos] = 0;
                    pos++;
                    left_index--;
                } else {
                    cout << "ERROR ON BUILD SPINE" << endl;
                }
            }
        }
    }
    //make dummy
    spine[pos]=1;
    pos++;
    spine[pos]=0;
    pos++;

#ifdef DEBUG
    assert(pos==spine_size+2);
#endif
}

template <typename T>double DFUDSEncoding::size(T dataStructure) {
    if(SIZE_MODE == 0)//get bpe
        return size_in_bytes(dataStructure)*8/(double)n;
    else if(SIZE_MODE == 1)
        return size_in_bytes(dataStructure)/double(1024);
    else if(SIZE_MODE == 2)
        return size_in_mega_bytes(dataStructure);
    else if(SIZE_MODE == 3)
        return size_in_mega_bytes(dataStructure)/double(1024);
}
double DFUDSEncoding::getSpaceOnlyBitVector() {
    double dfuds_size = dfuds.size();
    double spine_size = spine.size();
    double sum = dfuds_size+spine_size;
    return sum;
}
double DFUDSEncoding::getSpace() {
    double dfuds_size = size(dfuds);
    double spine_size = size(spine);
    double bps_size = size(bps);
    double depth_rrr_size = size(depthPos);
    double depth_element_size = depthValuesContainer.getSize()/double(n);
    double ldepth_element_size = lDepthValuesContainer.getSize()/double(n);
    double spine_rank_size = size(spineRank);
    double spine_select0_size = size(spine_ss_0);
    double spine_select1_size = size(spine_ss_1);
    double sum = dfuds_size+spine_size+bps_size+depth_element_size+depth_rrr_size+ldepth_element_size+spine_rank_size+spine_select0_size+spine_select1_size;
    return sum;
}
pair<DFUDSEncoding::arraysize,DFUDSEncoding::arraysize> DFUDSEncoding::top2_pos(DFUDSEncoding::arraysize i, DFUDSEncoding::arraysize j) {
    DFUDSEncoding::arraysize preorderu = rmq(i,j);
    DFUDSEncoding::arraysize r2m = 0;
    if(i==preorderu){
        r2m = rmq(preorderu+1,j);
    }else if(j==preorderu){
        r2m = rmq(i,preorderu-1);
    }else {
        DFUDSEncoding::arraysize u = arrayIndexToTreeIndex(preorderu + 1);
        int depthu = depth(u);
        int Ldepth = getLdepth(u);
        int L = degree(1);
        int lu = leftInnerSpineLengthNew(u, depthu);
        int Rdepth = depthu - 1;
        int Lleaves = preorderu - lLeavesRank.rank(u);
        int64_t spinePosition = 2 * preorderu - L - lu + Ldepth - Rdepth + 1 - Lleaves;
        DFUDSEncoding::arraysize v = arrayIndexToTreeIndex(rmq(i, preorderu - 1) + 1);
        DFUDSEncoding::arraysize w = arrayIndexToTreeIndex(rmq(preorderu + 1, j) + 1);
        int depthv = depth(v);
        int64_t lstart = depthv - depthu + spinePosition - spineRank.rank(spinePosition);
        int64_t rstart = degree(parent(w)) - childrank(w) - 1 + spineRank.rank(spinePosition);
        if (spine_ss_0.select(lstart + 1) < spine_ss_1.select(rstart + 1))
            r2m = rmq(i,preorderu-1);
        else
            r2m = rmq(preorderu+1,j);
    }
    return make_pair(preorderu,r2m);
}
pair<DFUDSEncoding::arraysize,DFUDSEncoding::arraysize> DFUDSEncoding::top2_pos_notsave(DFUDSEncoding::arraysize i,DFUDSEncoding::arraysize j)
{
    DFUDSEncoding::arraysize preorderu = rmq(i,j);
    DFUDSEncoding::arraysize r2m = 0;
    if(i==preorderu){
        r2m = rmq(preorderu+1,j);
    }else if(j==preorderu){
        r2m = rmq(i,preorderu-1);
    }else {
        DFUDSEncoding::arraysize u = arrayIndexToTreeIndex(preorderu + 1);
        int depthu = depthDecSequence(u);
        int Ldepth = lDepthDecSequence(u);
        int L = degree(1);
        int lu = leftInnerSpineLength(u);
        int Rdepth = depthu - 1;
        int Lleaves = preorderu - lLeavesRank.rank(u);
        int64_t spinePosition = 2 * preorderu - L - lu + Ldepth - Rdepth + 1 - Lleaves;
        DFUDSEncoding::arraysize v = arrayIndexToTreeIndex(rmq(i, preorderu - 1) + 1);
        DFUDSEncoding::arraysize w = arrayIndexToTreeIndex(rmq(preorderu + 1, j) + 1);
        int depthv = depthDecSequence(v);
        int64_t lstart = depthv - depthu + spinePosition - spineRank.rank(spinePosition);
        int64_t rstart = degree(parent(w)) - childrank(w) - 1 + spineRank.rank(spinePosition);
        if (spine_ss_0.select(lstart + 1) < spine_ss_1.select(rstart + 1))
            r2m = rmq(i,preorderu-1);
        else
            r2m = rmq(preorderu+1,j);
    }
    return make_pair(preorderu,r2m);
}
pair<DFUDSEncoding::arraysize,DFUDSEncoding::arraysize> DFUDSEncoding::top2_pos_rev_save(arraysize i, arraysize j)
{
    DFUDSEncoding::arraysize preorderu = rmq(i,j);
    DFUDSEncoding::arraysize r2m = 0;
    if(i==preorderu){
        r2m = rmq(preorderu+1,j);
    }else if(j==preorderu){
        r2m = rmq(i,preorderu-1);
    }else {
        DFUDSEncoding::arraysize u = arrayIndexToTreeIndex(preorderu + 1);
        int depthu = depth(u);
        int Ldepth = getLdepth(u);
        int lu = leftInnerSpineLengthNew(u, depthu);
        int L = degree(1);
        int Rdepth = depthu - 1;
        int Lleaves = preorderu - lLeavesRank.rank(u);
        int64_t spinePosition = 2 * preorderu - L - lu + Ldepth - Rdepth + 1 - Lleaves;
        DFUDSEncoding::arraysize v = arrayIndexToTreeIndex(rmq(i, preorderu - 1) + 1);
        DFUDSEncoding::arraysize w = arrayIndexToTreeIndex(rmq(preorderu + 1, j) + 1);
        int depthv = depth(v);
        int64_t lstart = depthv - depthu + spinePosition - spineRank.rank(spinePosition);
        int64_t rstart = degree(parent(w)) - childrank(w) - 1 + spineRank.rank(spinePosition);
        if (spine_ss_0.select(lstart + 1) <= spine_ss_1.select(rstart + 1))
            r2m = rmq(i,preorderu-1);
        else
            r2m = rmq(preorderu+1,j);
    }
    return make_pair(preorderu,r2m);
}
pair<DFUDSEncoding::arraysize,DFUDSEncoding::arraysize> DFUDSEncoding::top2_pos_rev_notsave(arraysize i, arraysize j)
{
    DFUDSEncoding::arraysize preorderu = rmq(i,j);
    DFUDSEncoding::arraysize r2m = 0;
    if(i==preorderu){
        r2m = rmq(preorderu+1,j);
    }else if(j==preorderu){
        r2m = rmq(i,preorderu-1);
    }else {
        DFUDSEncoding::arraysize u = arrayIndexToTreeIndex(preorderu + 1);
        int depthu = depthDecSequence(u);
        int Ldepth = lDepthDecSequence(u);
        int lu = leftInnerSpineLength(u);
        int L = degree(1);
        int Rdepth = depthu - 1;
        int Lleaves = preorderu - lLeavesRank.rank(u);
        int64_t spinePosition = 2 * preorderu - L - lu + Ldepth - Rdepth + 1 - Lleaves;
        DFUDSEncoding::arraysize v = arrayIndexToTreeIndex(rmq(i, preorderu - 1) + 1);
        DFUDSEncoding::arraysize w = arrayIndexToTreeIndex(rmq(preorderu + 1, j) + 1);
        int depthv = depthDecSequence(v);
        int64_t lstart = depthv - depthu + spinePosition - spineRank.rank(spinePosition);
        int64_t rstart = degree(parent(w)) - childrank(w) - 1 + spineRank.rank(spinePosition);
        if (spine_ss_0.select(lstart + 1) <= spine_ss_1.select(rstart + 1))
            r2m = rmq(i,preorderu-1);
        else
            r2m = rmq(preorderu+1,j);
    }
    return make_pair(preorderu,r2m);
}


DFUDSEncoding::arraysize DFUDSEncoding::rmqe(DFUDSEncoding::arraysize i, DFUDSEncoding::arraysize j) {//i<=x<=j
    if(i==j)
        return select(i+1);
    DFUDSEncoding::arraysize l = select(i+2);
    DFUDSEncoding::arraysize r = select(j+1);
    return bps.rmq_left(l,r);
}
DFUDSEncoding::arraysize DFUDSEncoding::rmq(DFUDSEncoding::arraysize i, DFUDSEncoding::arraysize j) {
    if(i==j)
        return i;
    DFUDSEncoding::arraysize w = rmqe(i,j);
    if(bps.rank0(bps.find_open(w)) == i+1){
        return i;
    }
    return bps.rank0(w)-1;
}
int DFUDSEncoding::getLdepth(DFUDSEncoding::arraysize x) {
#ifdef DEBUG
    assert(x>1);
#endif
    DFUDSEncoding::arraysize result = 0;
    DFUDSEncoding::arraysize curpos = x;
    while (curpos != 1 and depthPos[curpos] == 0) {
        result += degree(parent(curpos)) - childrank(curpos) - 1;
        curpos = parent(curpos);
    }
    if(curpos==1)
        return result;
    else
        //return result + lDepthValues[depthPosRank.rank(curpos)];
        return result + lDepthValuesContainer.get(depthPosRank.rank(curpos));
}
int DFUDSEncoding::getLdepth_init(DFUDSEncoding::arraysize x) {
#ifdef DEBUG
    assert(x>1 and x<bps.size());
#endif
    DFUDSEncoding::arraysize result = 0;
    DFUDSEncoding::arraysize curpos = x;
    while(curpos != 1){
        result += degree(parent(curpos)) - childrank(curpos) - 1;
        curpos = parent(curpos);
        if (depthPos[curpos] == 1) {
            result += lDepthValues[depthPosRank.rank(curpos)];
            break;
        }
    }
    return result;
}
DFUDSEncoding::arraysize DFUDSEncoding::leftInnerSpineLengthNew(DFUDSEncoding::arraysize u,int uDepth) {
#ifdef DEBUG
    assert(u>1);
#endif
    DFUDSEncoding::arraysize arrayIndexU = treeIndextoArrayIndex(u);
#ifdef DEBUG
    assert(arrayIndexU > 0);
#endif
    DFUDSEncoding::arraysize left_sibling = leftSibling(u);
    DFUDSEncoding::arraysize uMinusOne = arrayIndexToTreeIndex(arrayIndexU-1);
    int depthResult = 0;
    for(DFUDSEncoding::arraysize idx = uMinusOne;idx>left_sibling;idx=parent(idx)){
        if(depthPos[idx]==1){
            depthResult +=depthValues[depthPosRank.rank(idx)];
            return depthResult - uDepth + 1;
        }
        depthResult++;
    }
    return depthResult + 1;
}
DFUDSEncoding::arraysize DFUDSEncoding::leftInnerSpineLength(DFUDSEncoding::arraysize u) {
#ifdef DEBUG
    assert(u>1);
#endif
    int result = 0;
    if(childrank(u) == 0)
        return 0;
    DFUDSEncoding::arraysize left_sibling = leftSibling(u);
    for(DFUDSEncoding::arraysize i = left_sibling;not isLeaf(i);i = child(i,degree(i)-1)){
        result++;
    }
    return result+1;
}
DFUDSEncoding::arraysize DFUDSEncoding::depthDecSequence(DFUDSEncoding::arraysize x) {
    int result = 0;
    DFUDSEncoding::arraysize idx = x;
    while(idx!=1){
        idx = parent(idx);
        result++;
    }
    return result;
}
DFUDSEncoding::arraysize DFUDSEncoding::lDepthDecSequence(DFUDSEncoding::arraysize x) {
    DFUDSEncoding::arraysize result = 0;
    DFUDSEncoding::arraysize curpos = x;
    while (curpos != 1) {
        result += degree(parent(curpos)) - childrank(curpos) - 1;
        curpos = parent(curpos);
    }
    return result;
}
DFUDSEncoding::arraysize DFUDSEncoding::depth(DFUDSEncoding::arraysize x) {
    uint result = 0;
    DFUDSEncoding::arraysize idx = x;
    while(idx!=1){
        idx = parent(idx);
        result++;
        if(depthPos[idx]==1){
            result+=depthValuesContainer.get(depthPosRank.rank(idx));
            //result+=depthValues[depthPosRank.rank(idx)];
            break;
        }
    }
    return result;
}
DFUDSEncoding::arraysize DFUDSEncoding::depth_init(DFUDSEncoding::arraysize x) {//start at 0
#ifdef DEBUG
    assert(x<bps.size());
#endif
    int result = 0;
    DFUDSEncoding::arraysize idx = x;
    while(idx!=1){
        idx = parent(idx);
        result++;
        if(isCompleteDepth[idx]==1){
            result += allDepthValues[treeIndextoArrayIndex(idx)];
            break;
        }
    }
    allDepthValues.push_back(result);
    isCompleteDepth[x] = 1;
    return result;
}
DFUDSEncoding::arraysize DFUDSEncoding::leftSibling(DFUDSEncoding::arraysize x) {
#ifdef DEBUG
    assert(childrank(x)!=0 and x > 1);
#endif
    arraysize childrank_x = childrank(x);
    return child(parent(x),childrank_x-1);
}
DFUDSEncoding::arraysize DFUDSEncoding::parent(DFUDSEncoding::arraysize x) {
#ifdef DEBUG
    assert(x>1);
#endif
    //return max(prev(bps.find_open(x-1)),(arraysize)0)+1;
    return prev(bps.find_open(x-1))+1;
}
DFUDSEncoding::arraysize DFUDSEncoding::child(DFUDSEncoding::arraysize x,int i) {
#ifdef DEBUG
    assert(dfuds[x]==1);
#endif
    return bps.find_close(next(x)-(i+1))+1;
}
DFUDSEncoding::arraysize DFUDSEncoding::childrank(DFUDSEncoding::arraysize x) {
#ifdef DEBUG
    assert(x>1 and x<bps.size());
#endif
    DFUDSEncoding::arraysize y = bps.find_open(x-1);
    return next(y) - y - 1;
}
DFUDSEncoding::arraysize DFUDSEncoding::degree(DFUDSEncoding::arraysize x) {
#ifdef DEBUG
    assert(x<bps.size());
#endif
    if(isLeaf(x))
        return 0;
    return next(x) - x;
}
DFUDSEncoding::arraysize DFUDSEncoding::prev(DFUDSEncoding::arraysize x) {
#ifdef DEBUG
    assert(x<bps.size());
#endif
    DFUDSEncoding::arraysize rank0Value = bps.rank0(x);
    if(rank0Value <= 0)
        return 0;
    return select(rank0Value);
}
DFUDSEncoding::arraysize DFUDSEncoding::next(DFUDSEncoding::arraysize x) {
#ifdef DEBUG
    assert(x<bps.size());
#endif
    return select(bps.rank0(x)+1);
}
DFUDSEncoding::arraysize DFUDSEncoding::arrayIndexToTreeIndex(DFUDSEncoding::arraysize x) {
#ifdef DEBUG
    assert(x<bps.size()/2);
#endif
    if(x==0)
        return 1;
    return select(x)+1;
}
DFUDSEncoding::arraysize DFUDSEncoding::treeIndextoArrayIndex(DFUDSEncoding::arraysize x) {
#ifdef DEBUG
    assert(x<bps.size());
#endif
    return bps.rank0(x-1);
}
bool DFUDSEncoding::isLeaf(DFUDSEncoding::arraysize x) {
#ifdef DEBUG
    assert(x<bps.size());
#endif
    if(dfuds[x]==0)
        return true;
    else
        return false;
}
DFUDSEncoding::arraysize DFUDSEncoding::select(DFUDSEncoding::arraysize x){
    return bps.select(x);
}