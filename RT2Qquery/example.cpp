#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <stdlib.h>
#include "DFUDSEncoding.hpp"
#include "FischerHeunIndexing.hpp"
#include "DFUDSIndexing.hpp"
using namespace std;
using namespace sdsl;
typedef int arraysize;
void ourTest(string path,double usingSpace,int querySize);
void ourTestRev(string path,double usingSpace,int querySize);
template<typename T>
void Top2FromRMQTest(string path,int querySize);
void correctnessCheck(vector<arraysize>& A,arraysize start_index,arraysize end_index,arraysize rmq,arraysize r2m);
pair<arraysize,arraysize> second_largest_position_rmq_first(vector<arraysize>& arr,arraysize start,arraysize end);
vector<arraysize> readArrayFromFile(string PATH);
void printStat(string path,int querySize);
int ITERNUM = 1000;
int main(int argc,char** argv)
{
    srand(100);
    if(stoi(argv[1]) == 0){//Our RMQ TEST MODE
        if (argc < 5) {
            cout << "USAGE : ./test MODE PATH QUERYSIZE USINGSPACE" << endl;
            exit(0);
        }
        string path = argv[2];
        int querySize = stoi(argv[3]);
        double usingSpace = stod(argv[4]);
        ourTest(path,usingSpace,querySize);
    }else if(stoi(argv[1])==1) { // REVERSAL MODE
        if (argc < 6) {
            cout << "USAGE : ./test MODE PATH QUERYSIZE USINGSPACE" << endl;
            exit(0);
        }
        string path = argv[2];
        int querySize = stoi(argv[3]);
        double usingSpace = stod(argv[4]);
        ourTestRev(path,usingSpace,querySize);
    }else if(stoi(argv[1])==3){//TRIVIAL FISCHER TEST MODE
        if(argc < 4){
            cout << "USAGE : ./test MODE PATH QUERYSIZE" << endl;
            exit(0);
        }
        Top2FromRMQTest<FischerHeunIndexing>(argv[2],stoi(argv[3]));
    }else if(stoi(argv[1])==4){//STAT PRINT MODE
        if(argc < 4){
            cout << "USAGE : ./test MODE PATH QUERYSIZE" << endl;
            exit(0);
        }
        printStat(argv[2],stoi(argv[3]));
    }
}

void printStat(string path,int querySize){
    vector<int> A = readArrayFromFile(path);
    ofstream depthcsv;
    ofstream rmqdepthcsv;
    depthcsv.open("depth.csv");
    rmqdepthcsv.open("rmq"+to_string(querySize)+"depth.csv");
    GetInterval getInterval(A, 8,8,querySize);
    for(int i=1;i<getInterval.depth.size()+1;i++) {
        depthcsv << i << ",";
        rmqdepthcsv<< i <<",";
    }
    depthcsv<<"\n";
    rmqdepthcsv<<"\n";
    for(int i:getInterval.depth)
        depthcsv<<i<<",";
    for(int i:getInterval.rmqDepth)
        rmqdepthcsv<<i*100<<",";
    cout<<"DEPTH : "<<getInterval.depth.size()<<endl;
}
template<typename T>
void Top2FromRMQTest(string path,int querySize){
    vector<int> A = readArrayFromFile(path);
    T rmqStructure(A);
    cout<<"total space(bpe) : "<<rmqStructure.getSpaceRMQ()<<endl;
    double time = 0;
    for(int i=0;i<ITERNUM;i++) {
        int startIndex = rand() % (A.size() - querySize - 1);
        int endIndex = startIndex + querySize - 1;
        arraysize r2m;
        auto start_query = chrono::high_resolution_clock::now();
        arraysize rmqPos = rmqStructure.rmq(startIndex,endIndex);
        if(rmqPos == startIndex){
            r2m = rmqStructure.rmq(rmqPos+1,endIndex);
        }else if(rmqPos == endIndex){
            r2m = rmqStructure.rmq(startIndex,rmqPos-1);
        }else{
            arraysize left = rmqStructure.rmq(startIndex,rmqPos-1);
            arraysize right = rmqStructure.rmq(rmqPos+1,endIndex);
            r2m = (A[left] < A[right] ? right : left);
        }
        auto end_query = chrono::high_resolution_clock::now();
        correctnessCheck(A,startIndex,endIndex,rmqPos,r2m);
        time+= std::chrono::duration_cast<std::chrono::nanoseconds>(end_query - start_query).count() /(double) 1000;
        rmqStructure.initializeCache();
    }
    cout<<"elapsed time(ms) : "<<time/ITERNUM<<endl;
}
void ourTestRev(string path,double usingSpace,int querySize){
    vector<arraysize> A = readArrayFromFile(path);
    reverse(A.begin(),A.end());
    DFUDSEncoding* r2mq_rev;
    r2mq_rev = new DFUDSEncoding(A,usingSpace,true);
    reverse(A.begin(),A.end());
    cout<<"backward Mode"<<endl;
    cout<<"Query size : "<<querySize<<endl;
    cout<<"Using Space : "<<usingSpace<<endl;
    cout<<"total space(bpe) : "<<r2mq_rev->getSpace()<<endl;
    double time=0;
    for(int i=0;i<ITERNUM;i++) {
        int startIndex = rand() % (A.size() - querySize - 1);
        int endIndex = startIndex + querySize - 1;
        if(r2mq_rev->isSavingDepth) {
            auto start_query = chrono::high_resolution_clock::now();
            auto r2mqResult = r2mq_rev->top2_pos_rev_save(startIndex, endIndex);
            auto end_query = chrono::high_resolution_clock::now();
            time+= std::chrono::duration_cast<std::chrono::nanoseconds>(end_query - start_query).count() /(double) 1000;
            correctnessCheck(A, A.size()-endIndex-1, A.size()-startIndex-1,A.size()-1-r2mqResult.first, A.size()-1-r2mqResult.second);
        }else{
            auto start_query = chrono::high_resolution_clock::now();
            auto r2mqResult = r2mq_rev->top2_pos_rev_notsave(startIndex, endIndex);
            auto end_query = chrono::high_resolution_clock::now();
            time+= std::chrono::duration_cast<std::chrono::nanoseconds>(end_query - start_query).count() /(double) 1000;
            correctnessCheck(A, A.size()-endIndex-1, A.size()-startIndex-1,A.size()-1-r2mqResult.first, A.size()-1-r2mqResult.second);
        }
        r2mq_rev->bps.initializeCache();
    }
    cout<<"elapsed time(ms) : "<<time/ITERNUM<<endl;
    cout<<"--------------END EXPERIMENT----------------"<<endl;
}
void ourTest(string path,double usingSpace,int querySize){
    vector<arraysize> A = readArrayFromFile(path);
    DFUDSEncoding* r2mq;
    r2mq = new DFUDSEncoding(A,usingSpace,false);
    //cout<<"depth container size : "<<r2mq->depthValuesContainer.getSize()<<endl;
    cout<<"forward Mode"<<endl;
    cout<<"Query size : "<<querySize<<endl;
    cout<<"Using Space : "<<usingSpace<<endl;
    cout << "total space(bpe) : " << r2mq->getSpace() << endl;
    double time = 0;
    for (int i = 0; i < ITERNUM; i++) {
        int startIndex = rand() % (A.size() - querySize - 1);
        int endIndex = startIndex + querySize - 1;
        if(r2mq->isSavingDepth) {
            auto start_query = chrono::high_resolution_clock::now();
            auto r2mqResult = r2mq->top2_pos(startIndex, endIndex);
            auto end_query = chrono::high_resolution_clock::now();
            time += std::chrono::duration_cast<std::chrono::nanoseconds>(end_query - start_query).count() /(double) 1000;
            correctnessCheck(A, startIndex, endIndex, r2mqResult.first, r2mqResult.second);
        }else{
            auto start_query = chrono::high_resolution_clock::now();
            auto r2mqResult = r2mq->top2_pos_notsave(startIndex, endIndex);
            auto end_query = chrono::high_resolution_clock::now();
            time += std::chrono::duration_cast<std::chrono::nanoseconds>(end_query - start_query).count() /(double) 1000;
            correctnessCheck(A, startIndex, endIndex, r2mqResult.first, r2mqResult.second);
        }
        r2mq->bps.initializeCache();
    }
    cout << "elapsed time(ms) : " << time / ITERNUM << endl;
    cout<<"--------------END EXPERIMENT----------------"<<endl;
}
void correctnessCheck(vector<arraysize>& A,arraysize start_index,arraysize end_index,arraysize rmq,arraysize r2m)
{
    pair<arraysize,arraysize> real_top2 = second_largest_position_rmq_first(A,start_index,end_index);
    /*if(A[real_top2.first] != A[rmq])
        cout<<"RMQ ERROR ON : "<<start_index<<" , "<<end_index<<endl;
    if(A[real_top2.second] != A[r2m])
        cout<<"R2M ERROR ON : "<<start_index<<" , "<<end_index<<endl;*/
    if(real_top2.first != rmq)
        cout<<"RMQ ERROR ON : "<<start_index<<" , "<<end_index<<endl;
    if(real_top2.second != r2m)
        cout<<"R2M ERROR ON : "<<start_index<<" , "<<end_index<<endl;
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
pair<arraysize,arraysize> second_largest_position_rmq_first(vector<arraysize>& arr,arraysize start,arraysize end)//start <= target <= end, start starts 1
{
    arraysize i, first, second;
    arraysize first_pos, second_pos;
    if (arr.size() < 2)
    {
        cout<<" Invalid Input ";
        return make_pair(-1,-1);
    }
    first = second = 0;
    first_pos = second_pos = -1;
    for (i = start; i <= end ;i ++)
    {
        if (arr[i] > first)
        {
            second = first;
            first = arr[i];
            second_pos = first_pos;
            first_pos = i;
        }
        else if (arr[i] > second) {
            second = arr[i];
            second_pos = i;
        }
    }
    return make_pair(first_pos,second_pos);
}