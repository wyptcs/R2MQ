#include "FischerHeunIndexing.hpp"
#include <chrono>
#include <math.h>
#include <algorithm>
#include <vector>
#include <stack>
#include <stdint.h>
#include <cassert>
using namespace std;
typedef int16_t typeT;
typedef uint8_t typeP;
typedef int typeMDotDot;
typedef uint8_t typeMDot;
typeT FischerHeunIndexing::add0right(typeT x){
    return 2*x;
}
typeT FischerHeunIndexing::add1right(typeT x){
    return 2*x+1;
}
typeT FischerHeunIndexing::CartesianNumber(FischerHeunIndexing::arraysize i, FischerHeunIndexing::arraysize j, FischerHeunIndexing::arraysize *array) {
    typeT cartesian = 0;
    stack<FischerHeunIndexing::arraysize> stack;
    stack.push(i);
    cartesian = add1right(cartesian);
    for (int k = i+1; k <= j; k++) {
        while (!stack.empty() && array[k] > array[stack.top()]) {
            stack.pop();
            cartesian = add0right(cartesian);
        }
        stack.push(k);
        cartesian = add1right(cartesian);
    }
    while (stack.empty()) {
        stack.pop();
        cartesian = add0right(cartesian);
    }
    return cartesian;
}

void FischerHeunIndexing::initializeCartesians() {
    cartesians = new typeT[numOfBlocks];
    cartesianRMQs = new typeP[int(pow(2,2*sizeOfBlock))*sizeOfBlock*sizeOfBlock];
    for(int i=0;i<int(pow(2,2*sizeOfBlock))*sizeOfBlock*sizeOfBlock;i++)
        cartesianRMQs[i] = 0;
    for (int block = 0; block < numOfBlocks; block++) {
        int i = block*sizeOfBlock;
        int j = min(n-1, int((block+1)*sizeOfBlock - 1));
        typeT c = CartesianNumber(i, j,elements);
        cartesians[block] = c;
        create(cartesianRMQs,elements,c,i,j);
    }
}
void FischerHeunIndexing::create(typeP* RMQS,FischerHeunIndexing::arraysize* elements,FischerHeunIndexing::arraysize t,FischerHeunIndexing::arraysize i,FischerHeunIndexing::arraysize j){
    FischerHeunIndexing::arraysize index_t = t*sizeOfBlock*sizeOfBlock;
    FischerHeunIndexing::arraysize numElems = j - i + 1;
    for (FischerHeunIndexing::arraysize k = 0; k < numElems; k++) {
        RMQS[index_t+k*sizeOfBlock+k] = k;
    }
    for (FischerHeunIndexing::arraysize k = 0; k < numElems; k++) {
        for (FischerHeunIndexing::arraysize l = k+1; l < numElems; l++) {
            if (elements[RMQS[index_t+k*sizeOfBlock+l-1] + i] >= elements[i+l]) {
                RMQS[index_t+k*sizeOfBlock+l] = RMQS[index_t+k*sizeOfBlock+l-1];
            } else {
                RMQS[index_t+k*sizeOfBlock+l] = l;
            }
        }
    }
}
void FischerHeunIndexing::initializeCartesiansDot(){
    int tDotBlocks = int(ceil(double(numOfBlocks)/sizeOfBlock));
    cartesiansDot = new typeT[tDotBlocks];
    cartesianDotRMQs = new typeP[int(pow(2,2*sizeOfBlock))*sizeOfBlock*sizeOfBlock];
    for(int i=0;i<int(pow(2,2*sizeOfBlock))*sizeOfBlock*sizeOfBlock;i++)
        cartesianDotRMQs[i] = 0;
    for(int block = 0;block < tDotBlocks;block++){
        int i = block * sizeOfBlock;
        int j = min(numOfBlocks-1,(block+1)*sizeOfBlock-1);
        typeT c = CartesianNumber(i,j,ADot);
        cartesiansDot[block] = c;
        create(cartesianDotRMQs,ADot,c,i,j);
    }
}
FischerHeunIndexing::arraysize FischerHeunIndexing::notOrderedMinIndex(FischerHeunIndexing::arraysize index1,FischerHeunIndexing::arraysize index2) {
    FischerHeunIndexing::arraysize index;
    if(index1 <= index2)
        index = elements[index1] >= elements[index2] ? index1 : index2;
    else
        index = elements[index2] >= elements[index1] ? index2 : index1;
    return index;
}
FischerHeunIndexing::arraysize FischerHeunIndexing::minIndex(FischerHeunIndexing::arraysize index1, FischerHeunIndexing::arraysize index2) {
    return elements[index1] >= elements[index2] ? index1 : index2;
}
FischerHeunIndexing::arraysize FischerHeunIndexing::ADotIndextoAIndex(FischerHeunIndexing::arraysize i){
    FischerHeunIndexing::arraysize index;
    index = i*sizeOfBlock + cartesianRMQs[cartesians[i]*sizeOfBlock*sizeOfBlock+sizeOfBlock-1];
    return index;
}
FischerHeunIndexing::arraysize FischerHeunIndexing::minADotIndex(FischerHeunIndexing::arraysize index1,FischerHeunIndexing::arraysize index2){
    return ADot[index1] >= ADot[index2] ? index1 : index2;
}
void FischerHeunIndexing::initializeADot(){
    ADot = new FischerHeunIndexing::arraysize[numOfBlocks];
    for (int block = 0; block < numOfBlocks; block++) {
        FischerHeunIndexing::arraysize start = block*sizeOfBlock;
        FischerHeunIndexing::arraysize minidx = start;
        for (int i = start; i < min(int(start + sizeOfBlock), n); i++) {
            minidx = minIndex(minidx, i);
        }
        ADot[block] = elements[minidx];
    }
}
FischerHeunIndexing::arraysize FischerHeunIndexing::getMinIndexFromADot(FischerHeunIndexing::arraysize index1,FischerHeunIndexing::arraysize index2){//index1 <= i <= index2
    FischerHeunIndexing::arraysize result = index1;
    for(FischerHeunIndexing::arraysize i = index1+1;i <= index2;i++){
        result = minADotIndex(result,i);
    }
    return result;
}
void FischerHeunIndexing::buildMDotDot(){
    int nDot = numOfBlocks;
    int sDot = sizeOfBlock*sizeOfBlock;
    int rowSize = ceil(double(nDot)/double(sDot));
    int colSize = floor(log2(rowSize))+1;
    MDotDot = vector<vector<typeMDotDot>>(rowSize,vector<typeMDotDot>(colSize,0));
    for(int i = 0;i<rowSize;i++){
        MDotDot[i][0] = getMinIndexFromADot(i*sDot,min((i+1)*sDot-1,nDot-1));
    }
    for(int j=1;j<colSize;j++){
        int intervalLen = int(pow(2,j));
        for(int i = 0;i < rowSize - intervalLen + 1;i++){
            int startNext = i + intervalLen/2;
            MDotDot[i][j] = minADotIndex(MDotDot[i][j-1],MDotDot[startNext][j-1]);
        }
    }
}
void FischerHeunIndexing::buildMDot(){
    int nDot = numOfBlocks;
    int sDot = sizeOfBlock*sizeOfBlock;
    int s = sizeOfBlock;
    int rowSize = ceil(double(nDot)/double(s));
    int colSize = floor(log2(s))+1;
    vector<vector<int>> MDotTemp = vector<vector<int>>(rowSize,vector<int>(colSize,0));
    MDot = new typeMDot*[rowSize];
    for(int i=0;i<rowSize;i++)
        MDot[i] = new typeMDot[colSize];
    for(int i=0;i<rowSize;i++){
        MDotTemp[i][0] = getMinIndexFromADot(i*s,min((i+1)*s-1,nDot-1));
    }
    for(int j=1;j<colSize;j++){
        int intervalLen = int(pow(2,j));
        for(int i = 0;i < rowSize - intervalLen + 1;i++){
            int startNext = i + intervalLen/2;
            MDotTemp[i][j] = minADotIndex(MDotTemp[i][j-1],MDotTemp[startNext][j-1]);
        }
    }
    for(int i=0;i<rowSize;i++){
        for(int j=0;j<colSize;j++){
            MDot[i][j] = MDotTemp[i][j]-i*sizeOfBlock;
        }
    }
}
FischerHeunIndexing::arraysize FischerHeunIndexing::getOriginalMDotValue(FischerHeunIndexing::arraysize i,FischerHeunIndexing::arraysize j){
    FischerHeunIndexing::arraysize result;
    result = i*sizeOfBlock+MDot[i][j];
    return result;
}
FischerHeunIndexing::arraysize FischerHeunIndexing::bottomMin(FischerHeunIndexing::arraysize i, FischerHeunIndexing::arraysize j) {
    int iBlock = (int)(i/sizeOfBlock);
    int jBlock = (int)(j/sizeOfBlock);
    FischerHeunIndexing::arraysize end;
    if (iBlock == jBlock) {
        end = j;
    } else {
        end = (iBlock + 1)*sizeOfBlock - 1;
    }
    FischerHeunIndexing::arraysize firstMin = cartesianRMQs[cartesians[iBlock]*sizeOfBlock*sizeOfBlock+(i%sizeOfBlock)*sizeOfBlock+(end%sizeOfBlock)] + iBlock*sizeOfBlock;
    FischerHeunIndexing::arraysize start;
    if (iBlock == jBlock) {
        start = i;
    } else {
        start = jBlock*sizeOfBlock;
    }
    FischerHeunIndexing::arraysize secondMin = cartesianRMQs[cartesians[jBlock]*sizeOfBlock*sizeOfBlock+(start%sizeOfBlock)*sizeOfBlock+j%sizeOfBlock] + jBlock*sizeOfBlock;
    return minIndex(firstMin, secondMin);
}
FischerHeunIndexing::arraysize FischerHeunIndexing::ADotBottomMin(FischerHeunIndexing::arraysize topi,FischerHeunIndexing::arraysize topj){
    int topiBlock = topi/sizeOfBlock;
    int topjBlock = topj/sizeOfBlock;
#ifdef DEBUG
    assert(topiBlock <= topjBlock);
#endif
    if(topiBlock==topjBlock){
        FischerHeunIndexing::arraysize positionInADot = cartesianDotRMQs[cartesiansDot[topiBlock]*sizeOfBlock*sizeOfBlock+sizeOfBlock*(topi%sizeOfBlock)+topj%sizeOfBlock] + topiBlock*sizeOfBlock;
        return ADotIndextoAIndex(positionInADot);
    }else{
        FischerHeunIndexing::arraysize firstend = (topiBlock+1)*sizeOfBlock-1;
        FischerHeunIndexing::arraysize firstMin = cartesianDotRMQs[cartesiansDot[topiBlock]*sizeOfBlock*sizeOfBlock + sizeOfBlock*(topi%sizeOfBlock) + firstend%sizeOfBlock] + topiBlock*sizeOfBlock;
        FischerHeunIndexing::arraysize secondstart = topjBlock*sizeOfBlock;
        FischerHeunIndexing::arraysize secondMin = cartesianDotRMQs[cartesiansDot[topjBlock]*sizeOfBlock*sizeOfBlock + sizeOfBlock*(secondstart%sizeOfBlock) + topj%sizeOfBlock] + topjBlock*sizeOfBlock;
        return minIndex(ADotIndextoAIndex(firstMin),ADotIndextoAIndex(secondMin));
    }
}
FischerHeunIndexing::arraysize FischerHeunIndexing::ADotMiddleMin(int topiBlock,int topjBlock){
    int topiSuperBlock = topiBlock/sizeOfBlock;
    int topjSuperBlock = topjBlock/sizeOfBlock;
#ifdef DEBUG
        assert(topiSuperBlock <= topjSuperBlock);
#endif
        if(topiSuperBlock == topjSuperBlock){
            return getMinIndexFromTopBlock(topiBlock,topjBlock);
        }else{
            int firstend = (topiSuperBlock+1)*sizeOfBlock-1;
            FischerHeunIndexing::arraysize firstMinIndex = getMinIndexFromTopBlock(topiBlock,firstend);
            int secondstart = topjSuperBlock*sizeOfBlock;
            FischerHeunIndexing::arraysize secondMinIndex = getMinIndexFromTopBlock(secondstart,topjBlock);
            return minIndex(firstMinIndex,secondMinIndex);
        }
}
FischerHeunIndexing::arraysize FischerHeunIndexing::getMinIndexFromTopBlock(int topiBlock,int topjBlock){
    int interval = topjBlock - topiBlock;
    if(interval==0){
        FischerHeunIndexing::arraysize min = cartesianDotRMQs[cartesiansDot[topiBlock]*sizeOfBlock*sizeOfBlock+sizeOfBlock-1] + topiBlock*sizeOfBlock;
        return ADotIndextoAIndex(min);
    }else{
        int k = int(floor(log2(interval)));
        int twotok = int(pow(2,k));
        FischerHeunIndexing::arraysize firstMDotValue = getOriginalMDotValue(topiBlock,k);
        FischerHeunIndexing::arraysize secondMDotValue = getOriginalMDotValue(topjBlock-twotok+1,k);
        return minIndex(ADotIndextoAIndex(firstMDotValue),ADotIndextoAIndex(secondMDotValue));
    }
}
FischerHeunIndexing::arraysize FischerHeunIndexing::ADotTopMin(int topiSuperBlock,int topjSuperBlock){
    int interval = topjSuperBlock - topiSuperBlock;
    if(interval == 0){
        FischerHeunIndexing::arraysize MDotDotValue = MDotDot[topiSuperBlock][0];
        return ADotIndextoAIndex(MDotDotValue);
    }
    int k = int(floor(log2(interval)));
    int twotok = int(pow(2,k));
    FischerHeunIndexing::arraysize firstMDotDotValue = MDotDot[topiSuperBlock][k];
    FischerHeunIndexing::arraysize secondMDotDotValue = MDotDot[topjSuperBlock-twotok+1][k];
    return minIndex(ADotIndextoAIndex(firstMDotDotValue),ADotIndextoAIndex(secondMDotDotValue));
}
FischerHeunIndexing::FischerHeunIndexing(vector<FischerHeunIndexing::arraysize>& elems){
    n = elems.size();
    if (n == 0) return;
    elements = new FischerHeunIndexing::arraysize[n];
    for(FischerHeunIndexing::arraysize i=0;i<n;i++)
        elements[i] = elems[i];
    sizeOfBlock = int(log2(n)/4);
    if (sizeOfBlock < 1) return;
    numOfBlocks = (int)ceil((double)(n)/sizeOfBlock);
    initializeADot();
    buildMDotDot();
    buildMDot();
    initializeCartesians();
    initializeCartesiansDot();
    delete ADot;
}
double FischerHeunIndexing::getSpaceRMQ(){
    int rowSizeMDot = ceil(double(numOfBlocks)/double(sizeOfBlock));
    int colSizeMDot = floor(log2(sizeOfBlock))+1;
    int sizeT = numOfBlocks*sizeof(cartesians[0]);
    int sizeTDot = int(ceil(double(numOfBlocks)/sizeOfBlock))*sizeof(cartesiansDot[0]);
    int sizeP = int(pow(2,2*sizeOfBlock))*sizeOfBlock*sizeOfBlock;
    int sizePDot = int(pow(2,2*sizeOfBlock))*sizeOfBlock*sizeOfBlock;
    int sizeMDot = rowSizeMDot*colSizeMDot*sizeof(typeMDot);
    int sizeMDotDot = MDotDot.size()*MDotDot[0].size()*sizeof(typeMDotDot);
    int sum = sizeT + sizeTDot + sizeP + sizePDot + sizeMDot + sizeMDotDot;
    int sumBit = 8*sum;
    return sizeof(int)*8+double(sumBit)/double(n);
}
FischerHeunIndexing::arraysize FischerHeunIndexing::rmq(FischerHeunIndexing::arraysize i, FischerHeunIndexing::arraysize j) {
    FischerHeunIndexing::arraysize firstLayerRMQValue = bottomMin(i, j);
    int secondLayerStart = (int) (i / sizeOfBlock) + 1;
    int secondLayerEnd = (int) (j / sizeOfBlock) - 1;
    if (secondLayerEnd < secondLayerStart) return firstLayerRMQValue;
    FischerHeunIndexing::arraysize secondLayerRMQValue = notOrderedMinIndex(firstLayerRMQValue, ADotBottomMin(secondLayerStart, secondLayerEnd));
    int thirdLayerStart = secondLayerStart / sizeOfBlock + 1;
    int thirdLayerEnd = secondLayerEnd / sizeOfBlock - 1;
    if (thirdLayerEnd < thirdLayerStart) return secondLayerRMQValue;
    FischerHeunIndexing::arraysize thirdLayerRMQValue = notOrderedMinIndex(secondLayerRMQValue, ADotMiddleMin(thirdLayerStart, thirdLayerEnd));
    int fourthLayerStart = secondLayerStart / (sizeOfBlock*sizeOfBlock) + 1;
    int fourthLayerEnd = secondLayerEnd / (sizeOfBlock*sizeOfBlock) - 1;
    if (fourthLayerEnd < fourthLayerStart) return thirdLayerRMQValue;
    FischerHeunIndexing::arraysize fourthRMQValue = notOrderedMinIndex(thirdLayerRMQValue, ADotTopMin(fourthLayerStart, fourthLayerEnd));//TODO
    return fourthRMQValue;
    }
void FischerHeunIndexing::initializeCache(){
}

