#ifndef DFUDSENCODING_HPP
#define DFUDSENCODING_HPP

#include "sdsl_modified/bp_support_sada.hpp"
#include "sdsl_modified/rank_support_v5.hpp"
#include <sdsl/rrr_vector.hpp>
#include <sdsl/rmq_succinct_sada.hpp>

#include "BPIndexing.hpp"
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdint.h>
using namespace std;
using namespace sdsl;
enum MODE{CHAR=8,SHORT=16,INT=32};
class Depth
{
    MODE mode;
    vector<unsigned char> chardepth;
    vector<unsigned short> shortdepth;
    vector<unsigned int> intdepth;
public:
    Depth(){};
    Depth(MODE mode){
        this->mode = mode;
    }
    void push(int depth){
        if(mode==CHAR)
            chardepth.push_back(depth);
        else if(mode==SHORT)
            shortdepth.push_back(depth);
        else
            intdepth.push_back(depth);
    }
    int get(int i){
        if(mode==CHAR)
            return chardepth[i];
        else if(mode==SHORT)
            return shortdepth[i];
        else
            return intdepth[i];
    }
    int getSize(){
        if(mode==CHAR)
            return chardepth.size()*8;
        else if(mode==SHORT)
            return shortdepth.size()*16;
        else
            return intdepth.size()*32;
    }
};
class GetInterval
{
public:
    BPIndexing* r2mq;
    vector<int> depth;
    vector<int> rmqDepth;
    int arraySize;
    int depthVarSize;
    int lDepthVarSize;
public:
    GetInterval(vector<int>& array,int depthVarSize,int lDepthVarSize,int querySize)
    {
        r2mq = new BPIndexing(array);
        //cout<<"BUILD BP FOR CALCULATE DEPTH/LDEPTH STRUCTURE INTERVAL"<<endl;
        arraySize = array.size();
        int maxDepth = -1;
        for(int i=0;i<r2mq->bps.size();i++)
        {
            if(r2mq->bp[i] == 1){
                if(r2mq->bps.excess(i)-1 > maxDepth)
                    maxDepth = r2mq->bps.excess(i)-1;
            }
        }
        vector<int> depthStat(maxDepth+1,0);
        for(int i=0;i<r2mq->bps.size();i++)
        {
            if(r2mq->bp[i] == 1){
                depthStat[r2mq->bps.excess(i)-1]++;
            }
        }
        depth = depthStat;
        makeRmqDepthStat(array,querySize,maxDepth);
        this->depthVarSize = depthVarSize;
        this->lDepthVarSize = lDepthVarSize;
    }
    void makeRmqDepthStat(vector<int>& A,int querySize,int maxDepth){
        vector<int> rmqDepthStat(maxDepth+1,0);
        int iternum = 1000000;
        for(int i=0;i<iternum;i++) {
            int startIndex = rand() % (A.size() - querySize - 1);
            int endIndex = startIndex + querySize - 1;
            int rmqDepth = r2mq->bps.excess(r2mq->bps.rmq(r2mq->bps.select(startIndex + 2) - 1, r2mq->bps.select(endIndex + 2))) - 1;
            rmqDepthStat[rmqDepth]++;
        }
        /*for(int i=0;i<maxDepth+1;i++){
            cout<<rmqDepthStat[i]<<" ";
        }
        cout<<endl;
        for(int i=0;i<maxDepth+1;i++){
            cout<<depth[i]<<" ";
        }
        cout<<endl;*/
        rmqDepth = rmqDepthStat;
    }
    vector<int> getSavingCandidateOne(double usingSpace){
        vector<int> result;
        double currentUsingSpace = double(arraySize) / double((1.0/usingSpace) * (depthVarSize + lDepthVarSize) * 2);
        for(int i=(depth.size()-1)/2;i>=0;i--){
            if(depth[i] < currentUsingSpace){
                currentUsingSpace -= depth[i];
                result.push_back(i);
            }
        }
        reverse(result.begin(),result.end());
        return result;
    }
    vector<int> getSavingCandidateTwo(double usingSpace){
        vector<int> result;
        double currentUsingSpace = double(arraySize) / double((1.0/usingSpace) * (depthVarSize + lDepthVarSize) * 2);
        for(int i=0;i<depth.size();i++){
            if(depth[i] < currentUsingSpace){
                currentUsingSpace -= depth[i];
                result.push_back(i);
            }else{
                break;
            }
        }
        return result;
    }
    vector<int> getSavingCandidateThree(double usingSpace){
        vector<int> result;
        double currentUsingSpace = double(arraySize) / double((1.0/usingSpace) * (depthVarSize + lDepthVarSize) * 2);
        for(int i=1;i<depth.size();i=i*2){
            if(depth[i] < currentUsingSpace){
                currentUsingSpace -= depth[i];
                result.push_back(i);
            }else{
                break;
            }
        }
        return result;
    }
    vector<int> getSavingCandidateFour(double usingSpace){
        vector<int> result;
        int temp;
        double currentUsingSpace = double(arraySize) / double((1.0/usingSpace) * (depthVarSize + lDepthVarSize) * 2);
        for(int i=(depth.size()-1)/2;i>=0;i--){
            if(depth[i] < currentUsingSpace){
                temp = i;
                currentUsingSpace -= depth[i];
                result.push_back(i);
                break;
            }
        }
        for(int i=1;temp-i>=0;i=i*2){
            if(depth[temp-i] < currentUsingSpace){
                currentUsingSpace -= depth[temp-i];
                result.push_back(temp-i);
                temp = temp - i;
            }
        }
        reverse(result.begin(),result.end());
        return result;
    }
    vector<int> getSaving(double usingSpace){
        vector<int> result;
        double currentUsingSpace = double(arraySize) / double((1.0/usingSpace) * (depthVarSize + lDepthVarSize));
        for(int i=(depth.size()-1)/2;i>=0;i--){
            if(depth[i] < currentUsingSpace){
                currentUsingSpace -= depth[i];
                result.push_back(i);
            }
        }
        reverse(result.begin(),result.end());
        return result;
    }
    pair<int,int> getStartPosAndInterval(double usingSpace){
        int resultBlockSize = -1;
        int resultCurPos = -2;
        const double startDepthConst = 20;
        const double intervalConst = 40;
        if(depth.size() > int(intervalConst/usingSpace)*int(intervalConst/usingSpace)){
            int startDepth = int(startDepthConst/usingSpace);
            int interval = int(intervalConst/usingSpace);
            return make_pair(startDepth,interval);
        }
        for(int blockSize = depth.size()/2;blockSize > 1;blockSize--){
            int tempPos = getMinSumPos(usingSpace,blockSize);
            if(tempPos == -1){
                break;
            }
            resultCurPos = tempPos;
            resultBlockSize = blockSize;
        }
        if(resultCurPos==-2){
            for(int i=0;i<depth.size();i++){
                if(depth[i] > double(arraySize) / double((1.0/usingSpace) * (depthVarSize + lDepthVarSize))){
                    resultCurPos = i-1;
                    break;
                }
            }
            resultBlockSize = depth.size();
        }else{
            resultCurPos = getStarting(resultBlockSize,usingSpace);
        }
        return make_pair(resultCurPos,resultBlockSize);
    }
    int getStarting(int blockSize,double usingSpace){
        int starting = 0;
        int maxscore = -1;
        for(int i=0;i<blockSize;i++){
            if(startingScore(i,blockSize) > maxscore && isSpaceBoundSatisfy(i,blockSize,usingSpace)){
                starting = i;
                maxscore = startingScore(i,blockSize);
            }
        }
        return starting;
    }
    int startingScore(int starting, int blockSize){//DEBUG
        int score = 0;
        for(int i=0;i<starting;i++){
            score+=rmqDepth[i]*i;
        }
        for(int i=starting;i<rmqDepth.size();i++){
            score+=rmqDepth[i]*((i-starting)%blockSize);
        }
        return score;
    }
    bool isSpaceBoundSatisfy(int starting,int blockSize,double usingSpace){
        int currentResult = 0;
        for(int i=starting;i<depth.size();i=i+blockSize){
            currentResult += depth[i];
        }
        if (currentResult < double(arraySize) / double((1.0/usingSpace) * (depthVarSize + lDepthVarSize)))
            return true;
        else
            return false;
    }
    int getMinSumPos(double usingSpace, int blockSize){
        int curpos = -1;
        for(int i=0;i<blockSize;i++){
        //for(int i=blockSize-1;i>=0;i--){
            int currentResult = 0;
            for(int j=i;j<depth.size();j=j+blockSize){
                currentResult += depth[j];
            }
            if (currentResult < double(arraySize) / double((1.0/usingSpace) * (depthVarSize + lDepthVarSize)))
                curpos = i;
        }
        return curpos;
    }
};
class DFUDSEncoding {
    typedef int arraysize;
    static const int BLOCK_SIZE = 512;
private:

    int n;
    bool reverse;


    vector<int> allDepthValues;//TEMPORARY USED
    bit_vector isCompleteDepth;//TEMPORARY USED
    vector<int> allLDepthValues;//TEMPORARY USED
    bit_vector isCompleteLDepth;//TEMPORARY USED

    rank_support_v5<100,3> lLeavesRank;
    rank_support_v5<> spineRank;
    select_support_mcl<0,1> spine_ss_0;
    select_support_mcl<1,1> spine_ss_1;

    Depth lDepthValuesContainer;

    vector<arraysize> depthValues;
    rrr_vector<> depthPos;
    rrr_vector<>::rank_1_type depthPosRank;
    vector<arraysize> lDepthValues;
    void makeRMQStructure(vector<arraysize>& array);
    void makeTop2Structure(vector<arraysize>& array);

    void makeDepthValuesContainer(MODE mode);
    void makeLDepthValuesContainer(MODE mode);
    void makeDepthPosNew(arraysize n,vector<int> savePoint);
    void makeLDepthPosNew(arraysize n,vector<int> savePoint);
    void makeDepthPos(arraysize n,int startDepth,int interval);
    void makeLDepthPos(arraysize n, int startDepth, int interval);
    arraysize rmqe(arraysize i,arraysize j);
    int getLdepth_init(arraysize x);
    int getLdepth(arraysize x);
    arraysize select(arraysize x);
    template <typename T> double size(T dataStructure);

    arraysize depthDecSequence(arraysize x);
    arraysize lDepthDecSequence(arraysize x);


    arraysize depth_init(arraysize x);
    arraysize leftSibling(arraysize x);
    arraysize parent(arraysize x);
    arraysize childrank(arraysize x);//childrank start at 0
    arraysize child(arraysize x,int i);//child index start at 0
    arraysize degree(arraysize x);
    arraysize prev(arraysize x);//when dfuds[x] == 0, return current pos
    arraysize next(arraysize x);//when dfuds[x] == 0, return next pos

    arraysize treeIndextoArrayIndex(arraysize x);
    bool isLeaf(arraysize x);

public:
    bit_vector dfuds;
    bit_vector spine;
    int maxDepth;
    int interval;
    bool isSavingDepth;
    Depth depthValuesContainer;
    static const int SIZE_MODE = 0;//0 : byte, 1 : KB, 2 : MB 3 : GB
    bp_support_sada<BLOCK_SIZE,1,rank_support_v5<>,select_support_mcl<0,1>> bps;
    DFUDSEncoding(vector<arraysize>&array,double usingSpace,bool reverse);
    DFUDSEncoding(vector<arraysize>&array);
    DFUDSEncoding();

    void makeDFUDS(vector<arraysize>& array);
    void makeSpineUsingTwoDVector(vector<arraysize>& array);
    pair<arraysize,arraysize> top2_pos(arraysize i,arraysize j);
    pair<arraysize,arraysize> top2_pos_notsave(arraysize i,arraysize j);
    pair<arraysize,arraysize> top2_pos_rev_save(arraysize i,arraysize j);
    pair<arraysize,arraysize> top2_pos_rev_notsave(arraysize i,arraysize j);
    arraysize rmq(arraysize i,arraysize j);
    arraysize leftInnerSpineLength(arraysize u);
    arraysize leftInnerSpineLengthNew(arraysize u,int uDepth);
    double getSpace();
    double getSpaceOnlyBitVector();
    double getSpaceHuffman();

    arraysize depth(arraysize x);//TEMP PUBLIC
    arraysize arrayIndexToTreeIndex(arraysize x);//TEMP PUBLIC
};

#endif
