#ifndef FISCHERHEUNINDEXING_HPP
#define FISCHERHEUNINDEXING_HPP
#include <chrono>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stack>
#include <stdint.h>
using namespace std;
typedef int16_t typeT;
typedef uint8_t typeP;
typedef int typeMDotDot;
typedef uint8_t typeMDot;
class FischerHeunIndexing{
    typedef int arraysize;
private:
    typeT* cartesians;//T
    typeT* cartesiansDot;//t'
    arraysize* elements;
    typeP* cartesianRMQs;
    typeP* cartesianDotRMQs;
    vector<vector<typeMDotDot>> MDotDot;
    typeMDot** MDot;
    arraysize n;  // size of array
    int sizeOfBlock;  // size of blocks
    int numOfBlocks;
    arraysize* ADot;
    typeT add0right(typeT x);
    typeT add1right(typeT x);
    typeT CartesianNumber(arraysize i,arraysize j,arraysize* array);
    void initializeCartesians();
    void create(typeP* RMQS,arraysize* elements,arraysize t,arraysize i,arraysize j);
    void initializeCartesiansDot();
    arraysize notOrderedMinIndex(arraysize index1,arraysize index2);
    arraysize minIndex(arraysize index1, arraysize index2);
    arraysize ADotIndextoAIndex(arraysize i);
    arraysize minADotIndex(arraysize index1,arraysize index2);
    void initializeADot();
    arraysize getMinIndexFromADot(arraysize index1,arraysize index2);
    void buildMDotDot();
    void buildMDot();
    arraysize getOriginalMDotValue(arraysize i,arraysize j);
    arraysize bottomMin(arraysize i, arraysize j);
    arraysize ADotBottomMin(arraysize topi,arraysize topj);
    arraysize ADotMiddleMin(int topiBlock,int topjBlock);
    arraysize getMinIndexFromTopBlock(int topiBlock,int topjBlock);
    arraysize ADotTopMin(int topiSuperBlock,int topjSuperBlock);
public:
    FischerHeunIndexing(vector<arraysize>& elems);
    double getSpaceRMQ();
    arraysize rmq(arraysize i, arraysize j);
    void initializeCache();
};
#endif
