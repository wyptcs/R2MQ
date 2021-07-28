#include <iostream>
#include <chrono>
#include <vector>
#include <string>
#include <stdlib.h>
//#include "sdsl_SEA2017/include/sdsl/rmq_succinct_rec.hpp"
#include "sdsl_SEA2017/include/sdsl/rmq_succinct_rec_new.hpp"
using namespace std;
using namespace sdsl;
typedef int arraysize;
void SEARMQTest(string path,int querySize);
void correctnessCheck(vector<arraysize>& A,arraysize start_index,arraysize end_index,arraysize rmq,arraysize r2m);
pair<arraysize,arraysize> second_smallest_position_rmq_first(vector<arraysize>& arr,arraysize start,arraysize end);
vector<arraysize> readArrayFromFile(string PATH);
void printStat(string path,int querySize);
int ITERNUM = 1000;
int main(int argc,char** argv) {
    srand(100);
    SEARMQTest(argv[1],stoi(argv[2]));
}
void SEARMQTest(string path,int querySize){
    vector<arraysize> A = readArrayFromFile(path);
    for(int i=0;i<A.size();i++)
        A[i] = -A[i];
    rmq_succinct_rec<1024,128,true> rmq(&A);
    cout<<"total space without array(bpe)"<<double(size_in_bytes(rmq))*8/double(A.size())<<endl;
    double time = 0;
    for(int i=0;i<ITERNUM;i++) {
        int startIndex = rand() % (A.size() - querySize - 1);
        int endIndex = startIndex + querySize - 1;
        arraysize r2m;
        auto start_query = chrono::high_resolution_clock::now();
        arraysize rmqPos = rmq(startIndex,endIndex);
        //auto end_query = chrono::high_resolution_clock::now();
        if(rmqPos == startIndex){
            r2m = rmq(rmqPos+1,endIndex);
        }else if(rmqPos == endIndex){
            r2m = rmq(startIndex,rmqPos-1);
        }else{
            arraysize left = rmq(startIndex,rmqPos-1);
            arraysize right = rmq(rmqPos+1,endIndex);
            r2m = (A[left] > A[right] ? right : left);
        }
        auto end_query = chrono::high_resolution_clock::now();
        correctnessCheck(A,startIndex,endIndex,rmqPos,r2m);
        time+= std::chrono::duration_cast<std::chrono::nanoseconds>(end_query - start_query).count() /(double) 1000;

    }
    cout<<"elapsed time(ms) : "<<time/ITERNUM<<endl;
}
void correctnessCheck(vector<arraysize>& A,arraysize start_index,arraysize end_index,arraysize rmq,arraysize r2m)
{
    pair<arraysize,arraysize> real_top2 = second_smallest_position_rmq_first(A,start_index,end_index);
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
pair<arraysize,arraysize> second_smallest_position_rmq_first(vector<arraysize>& arr,arraysize start,arraysize end)//start <= target <= end, start starts 1
{
    arraysize i, first, second;
    arraysize first_pos, second_pos;
    if (arr.size() < 2)
    {
        cout<<" Invalid Input ";
        return make_pair(-1,-1);
    }
    first = second = INT32_MAX;
    first_pos = second_pos = -1;
    for (i = start; i <= end ;i ++)
    {
        if (arr[i] < first)
        {
            second = first;
            first = arr[i];
            second_pos = first_pos;
            first_pos = i;
        }
        else if (arr[i] < second) {
            second = arr[i];
            second_pos = i;
        }
    }
    return make_pair(first_pos,second_pos);
}