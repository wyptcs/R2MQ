#include <iostream>
#include <vector>
#include <algorithm>
#include <bits/stdc++.h>
#include <stack>
#include <chrono>
#include <unordered_map>
#include <sdsl/rrr_vector.hpp>
#include <string>
using namespace std;
using namespace sdsl;
bit_vector tobitvector(vector<bool> input);
struct Node_huff
{
    char ch;
    int freq;
    Node_huff *left, *right;
};
struct comp
{
    bool operator()(Node_huff* l, Node_huff* r)
    {
        // highest priority item has lowest frequency
        return l->freq > r->freq;
    }
};
Node_huff* getNode(char ch, int freq, Node_huff* left, Node_huff* right);
void encode(Node_huff* root, string str,unordered_map<char, string> &huffmanCode);
void decode(Node_huff* root, int &index, string str);
string buildHuffmanTree(string text);
int getHuffmanSize(bit_vector* b);
int getLz77Size(bit_vector* b);
pair<int,int> second_largest_position_rmq_first(vector<int>& arr,int start,int end);
vector<pair<double,double>> dagTestUsingOptimization(vector<int>& example);
vector<pair<double,double>> dagTest(vector<int>& example);
vector<pair<double,double>> spineTest(vector<int>& example,bool reverse);