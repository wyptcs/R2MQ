#include <iostream>
#include <chrono>
#include <algorithm>
#include <cassert>
#include <vector>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iterator>
#include <sdsl/rrr_vector.hpp>
#include "topk_third.cpp"
using namespace std;
using namespace sdsl;
vector<int> readArrayFromFile(string PATH);
class Test
{
private:        
    Pawel* pawel;
public:
    Test()
    {

    }
    void printDagResult(vector<int>& array,bool isOptimization)
    {   
        vector<pair<double,double>> result;
        if(isOptimization){
	result = dagTestUsingOptimization(array);
        } else {
	result = dagTest(array);
        }
        cout<<"SPACE OF DAG STRUCTURE(bpe) : "<<result[0].first<<endl;        
        cout<<"DAG CONSTRUCTION TIME : "<<result[0].second<<endl;
        cout<<"DAG STRUCTURE COMPRESSED USING HUFFMAN(bpe) : "<<result[1].first<<endl;
        cout<<"HUFFMAN COMPRESSION TIME : "<<result[1].second<<endl;
        cout<<"DAG STRUCTURE COMPRESSED USING LZW(bpe) : "<<result[2].first<<endl;
        cout<<"LZW COMPRESSION TIME : "<<result[2].second<<endl;
    }
    void printSpineResult(vector<int>& array,bool revMode)
    {
        auto result = spineTest(array,revMode);
        cout<<"SPACE OF SPINE STRUCTURE(bpe) : "<<result[0].first<<endl;
        cout<<"SPINE CONSTRUCTION TIME : "<<result[0].second<<endl;
        cout<<"SPINE STRUCTURE COMPRESSED USING HUFFMAN(bpe) : "<<result[1].first<<endl;
        cout<<"HUFFMAN COMPRESSION TIME : "<<result[1].second<<endl;
        cout<<"SPINE STRUCTURE COMPRESSED USING LZW(bpe) : "<<result[2].first<<endl;
        cout<<"LZW COMPRESSION TIME : "<<result[2].second<<endl;
    }
    void printPawelResult(vector<int>& array)
    {
        pawel = new Pawel(array,2);
        auto result = pawel->getResult();
        cout<<"SPACE OF OPTIMAL STRUCTURE(bpe) : "<<result[0].first<<endl;        
        cout<<"OPTIMAL CONSTRUCTION TIME : "<<result[0].second<<endl;
        cout<<"OPTIMAL STRUCTURE COMPRESSED USING HUFFMAN(bpe) : "<<result[1].first<<endl;
        cout<<"HUFFMAN COMPRESSION TIME : "<<result[1].second<<endl;
        cout<<"OPTIMAL STRUCTURE COMPRESSED USING LZW(bpe) : "<<result[2].first<<endl;
        cout<<"LZW COMPRESSION TIME : "<<result[2].second<<endl;
    }
};
vector<int> readArrayFromFile(string PATH)
{
    ifstream inFile(PATH);
    vector<int> arr;
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
int main(int argc,char** argv)
{
    string path = argv[1];
    bool optMode = stoi(argv[2]) == 1 ? true : false;
    bool revMode = stoi(argv[3]) == 1 ? true : false;
    vector<int> A = readArrayFromFile(path);
    Test* t = new Test();
    t->printDagResult(A,optMode);
    cout<<endl;
    cout<<endl;
    t->printPawelResult(A);
    cout<<endl;
    cout<<endl;
    t->printSpineResult(A,revMode);
    return 0;
}