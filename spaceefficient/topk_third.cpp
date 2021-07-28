#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <math.h>
#include <utility>
#include <cmath>
#include <list>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <sdsl/wavelet_trees.hpp>
#include <queue>
#include <iostream>
#include <string>
#include <unordered_map>
#include <array>
#include <chrono>
#include "auxiliary.hpp"
using namespace std;
using namespace sdsl;
class Pawel{
private:
    rrr_vector<255> S;
    int huffmanSize;
    int lzwSize;
    int originalSize;
    double buildTime;
    double huffmanTime;
    double lzwTime;
public:
    vector<int> startingpointvector;
    vector<int> numberOfTwo;
    long int comparecount;
    Pawel(){

    }
    Pawel(const vector<int> input,int k){
        bit_vector SS;
        comparecount = 0;
        int ele = 0;
        int n = input.size();
        int delta = 0;
        int sigma = 0;
        int sigma_prev = 0;
        int pivot = 0;
        int startingpoint = 0;
        SS.resize((k+1)*n);
        auto start_encoding = chrono::high_resolution_clock::now();
        vector<int> S_prev(n,0);
        vector<int> S_curr(n,0);
        for(int i=1;i<n;i++){
            for(int j=startingpoint;j<i;j++){
                if(S_curr[j] < k && input[j] < input[i]){
                    S_curr[j]++;
                    sigma +=1;
                }
            }
            for(int j=startingpoint;j<i;j++){
                if(S_curr[j] < k)
                    break;
                else
                    startingpoint++;
            }
            int_to_unary(SS,sigma,&pivot);
            sigma=0;
        }
        auto end_encoding = chrono::high_resolution_clock::now();
        buildTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_encoding - start_encoding).count() /(double) 1000;
        SS.resize(pivot);
        originalSize = SS.size();
        huffmanSize = getHuffmanSize(&SS);
        auto end_huffman = chrono::high_resolution_clock::now();
        huffmanTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_huffman - end_encoding).count() /(double) 1000;
        lzwSize = getLz77Size(&SS);
        auto end_lzw = chrono::high_resolution_clock::now();
        lzwTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end_lzw - end_huffman).count() /(double) 1000;
        rrr_vector<255> rrr_S(SS);
        S=rrr_S;        
    }
    vector<pair<double,double>> getResult(){
        vector<pair<double,double>> result;
        result.push_back(make_pair(originalSize,buildTime));
        result.push_back(make_pair(huffmanSize,huffmanTime));
        result.push_back(make_pair(lzwSize,lzwTime));
        return result;
    }
    void int_to_unary(bit_vector& S,int input,int* pivot){
        for(int i=0;i<input;i++){
            S[*pivot] = 1;
            (*pivot)++;
        }
        S[*pivot] = 0;
        (*pivot)++;
    }
};