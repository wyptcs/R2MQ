#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <vector>
#include <bitset>
#include "sdsl_SEA2017/include/sdsl/rmq_succinct_rec_new.hpp"
#include <string>
using namespace std;
typedef uint arraysize;
vector<arraysize> readArrayFromFile(string PATH);
vector<bool> createIncreaseBitSequence(vector<arraysize>& arr);
vector<bool> createDecreaseBitSequence(vector<arraysize>& arr);
vector<bool> createIncDecBitSequence(vector<arraysize>& arr);
bool initMode(vector<arraysize>& arr,int i);
int main(int argc,char** argv)
{
    vector<arraysize> testarray = readArrayFromFile(argv[1]);
    auto result_incdec = createIncDecBitSequence(testarray);
    auto result_inc = createIncreaseBitSequence(testarray);
    auto result_dec = createDecreaseBitSequence(testarray);
    int elemnum_incdec = count(result_incdec.begin(),result_incdec.end(),true);
    int elemnum_inc = count(result_inc.begin(),result_inc.end(),true);
    int elemnum_dec = count(result_dec.begin(),result_dec.end(),true);
    cout<<"element number_inc : "<<double(elemnum_inc)*32/double(testarray.size())<<endl;
    cout<<"element number_dec : "<<double(elemnum_dec)*32/double(testarray.size())<<endl;
    cout<<"element number_incdec : "<<double(elemnum_incdec)*32/double(testarray.size())<<endl;

    /*for(int i=0;i<result_inc.size();i++){
        if(result_inc[i])
            cout<<1;
        else
            cout<<0;
    }
    cout<<endl;
    for(int i=0;i<result_dec.size();i++){
        if(result_dec[i])
            cout<<1;
        else
            cout<<0;
    }
    cout<<endl;
    for(int i=0;i<result_incdec.size();i++){
        if(result_incdec[i])
            cout<<1;
        else
            cout<<0;
    }
    cout<<endl;*/
    //vector<arraysize> toy = {7,8,9,10,4,3,2,1};
    //vector<arraysize> toy = {7,8,9,10,1,2,3,4,5,3};
    //vector<arraysize> toy = {4,5,6,3,2};
    //vector<arraysize> toy = {4,3,2,1,10,9,8,7};
    //auto result = createIncreaseBitSequence(toy);
    //auto result = createDecreaseBitSequence(toy);
    //auto result = createIncDecBitSequence(toy);
    /*for(int i=0;i<result.size();i++){
        if(result[i])
            cout<<1;
        else
            cout<<0;
    }
    cout<<endl;*/
}
bool initMode(vector<arraysize>& arr,int i){
    if(arr[i] < arr[i+1])
        return true;
    else
        return false;
}
vector<bool> createIncDecBitSequence(vector<arraysize>& arr){
    int n = arr.size();
    vector<bool> resultBitSequence(n,false);
    int pivot = 0;
    bool isIncrMode = initMode(arr,0);
    for(int i=0;i<n-1;i++){
        if(arr[i] < arr[i+1] and isIncrMode)
            continue;
        else if(arr[i] > arr[i+1] and (not isIncrMode))
            continue;
        else{
            if(isIncrMode){
                resultBitSequence[i] = true;
                resultBitSequence[i-1] = true;
            }else{
                resultBitSequence[pivot] = true;
                resultBitSequence[pivot+1] = true;
            }
            pivot = i+1;
            if(pivot < n-1){
                isIncrMode = initMode(arr,pivot);
            }
        }
    }
    if(isIncrMode){
        resultBitSequence[n-2] = true;
        resultBitSequence[n-1] = true;
    }else{
        resultBitSequence[pivot] = true;
        resultBitSequence[pivot+1] = true;
    }
    return resultBitSequence;
}
vector<bool> createDecreaseBitSequence(vector<arraysize>& arr){
    int n = arr.size();
    vector<bool> resultBitSequence(n,false);
    int pivot = 0;
    for(int i=0;i<n-1;i++){
        if(arr[i] > arr[i+1])
            continue;
        else{
            resultBitSequence[pivot] = true;
            if(pivot < i)
                resultBitSequence[pivot+1] = true;
            pivot = i+1;
        }
    }
    if(arr[n-2] > arr[n-1])
    {
        resultBitSequence[pivot] = true;
        resultBitSequence[pivot+1] = true;
    }else{
        resultBitSequence[n-1] = true;
    }
    return resultBitSequence;
}
vector<bool> createIncreaseBitSequence(vector<arraysize>& arr)
{
    int n = arr.size();
    vector<bool> resultBitSequence(n,false);
    int pivot = 0;
    for(int i=0;i<n-1;i++){
        if(arr[i] < arr[i+1])
            continue;
        else{
            resultBitSequence[i] = true;
            if(pivot < i)
                resultBitSequence[i-1] = true;
            pivot = i+1;
        }
    }
    if(arr[n-2] < arr[n-1])
        resultBitSequence[n-2] = true;
    resultBitSequence[n-1] = true;
    return resultBitSequence;
}
vector<arraysize> readArrayFromFile(string PATH)
{
    ifstream inFile(PATH);
    vector<arraysize> arr;
    if(inFile.is_open())
    {
        while(not inFile.eof())
        {
            string str;
            getline(inFile, str);
            if(str.size()!=0)
                arr.push_back(stoi(str));
        }
    }
    return arr;
}