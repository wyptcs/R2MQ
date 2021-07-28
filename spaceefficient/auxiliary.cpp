#include <iostream>
#include <vector>
#include <algorithm>
#include <bits/stdc++.h>
#include <stack>
#include <queue>
#include <set>
#include <chrono>
#include <string>
#include <unordered_map>
#include <sdsl/rrr_vector.hpp>
#include <iterator>
#include "newlz.cpp"
#include "sdsl_SEA2017/include/sdsl/rmq_succinct_rec_new.hpp"
using namespace std;
using namespace sdsl;
bit_vector tobitvector(vector<bool> input){
    bit_vector result;
    result.resize(input.size());
    for(int i=0;i<input.size();i++){
        if(input[i])
            result[i] = 1;
        else
            result[i] = 0;
    }
    return result;
}
struct Node_huff
{
    char ch;
    int freq;
    Node_huff *left, *right;
};

// Function to allocate a new tree node
Node_huff* getNode(char ch, int freq, Node_huff* left, Node_huff* right)
{
    Node_huff* node = new Node_huff();

    node->ch = ch;
    node->freq = freq;
    node->left = left;
    node->right = right;

    return node;
}
// Comparison object to be used to order the heap
struct comp
{
    bool operator()(Node_huff* l, Node_huff* r)
    {
        // highest priority item has lowest frequency
        return l->freq > r->freq;
    }
};

// traverse the Huffman Tree and store Huffman Codes
// in a map.
void encode(Node_huff* root, string str,
            unordered_map<char, string> &huffmanCode)
{
    if (root == nullptr)
        return;

    // found a leaf node
    if (!root->left && !root->right) {
        huffmanCode[root->ch] = str;
    }

    encode(root->left, str + "0", huffmanCode);
    encode(root->right, str + "1", huffmanCode);
}

// traverse the Huffman Tree and decode the encoded string
void decode(Node_huff* root, int &index, string str)
{
    if (root == nullptr) {
        return;
    }

    // found a leaf node
    if (!root->left && !root->right)
    {
        cout << root->ch;
        return;
    }

    index++;

    if (str[index] =='0')
        decode(root->left, index, str);
    else
        decode(root->right, index, str);
}

// Builds Huffman Tree and decode given input text
string buildHuffmanTree(string text)
{
    // count frequency of appearance of each character
    // and store it in a map

    unordered_map<char, int> freq;
    for (char ch: text) {
        freq[ch]++;
    }

    // Create a priority queue to store live nodes of
    // Huffman tree;
    priority_queue<Node_huff*, vector<Node_huff*>, comp> pq;

    // Create a leaf node for each character and add it
    // to the priority queue.
    for (auto pair: freq) {
        pq.push(getNode(pair.first, pair.second, nullptr, nullptr));
    }

    // do till there is more than one node in the queue
    while (pq.size() != 1)
    {
        // Remove the two nodes of highest priority
        // (lowest frequency) from the queue
        Node_huff *left = pq.top(); pq.pop();
        Node_huff *right = pq.top();	pq.pop();

        // Create a new internal node with these two nodes
        // as children and with frequency equal to the sum
        // of the two nodes' frequencies. Add the new node
        // to the priority queue.
        int sum = left->freq + right->freq;
        pq.push(getNode('\0', sum, left, right));
    }

    // root stores pointer to root of Huffman Tree
    Node_huff* root = pq.top();

    // traverse the Huffman Tree and store Huffman Codes
    // in a map. Also prints them
    unordered_map<char, string> huffmanCode;
    encode(root, "", huffmanCode);

    // print encoded string
    string str = "";
    for (char ch: text) {
        str += huffmanCode[ch];
    }
    return str;
}
int getHuffmanSize(bit_vector* b){
    vector<char> S_new;
    int block_num = ((*b).size() - 1) / 8 + 1;
    (*b).resize(block_num*8);
    for(int i = (*b).size();i<block_num*8;i++)
        (*b)[i] = 1;//DUMMY
    for (int i = 0; i < block_num; i++) {
        char result_char = 0;
        for (int j = 7; j >= 0; j--) {
            result_char += ((*b)[i * 8 + 7 - j] << j);
        }
        S_new.push_back(result_char);
    }
    string st(S_new.begin(), S_new.end());
    string s_huff = buildHuffmanTree(st);
    return s_huff.size();
}
int getLz77Size(bit_vector* b){
    vector<char> S_new;
    int block_num = ((*b).size() - 1) / 8 + 1;
    (*b).resize(block_num*8);
    for(int i = (*b).size();i<block_num*8;i++)
        (*b)[i] = 1;//DUMMY
    for (int i = 0; i < block_num; i++) {
        char result_char = 0;
        for (int j = 7; j >= 0; j--) {
            result_char += ((*b)[i * 8 + 7 - j] << j);
        }
        S_new.push_back(result_char);
    }
    string st(S_new.begin(), S_new.end());
    vector<int> compressed = compress(st);
    return compressed.size()*8;
}
pair<int,int> second_largest_position_rmq_first(vector<int>& arr,int start,int end)//start <= target <= end, start starts 1
{
    int i, first, second;
    int first_pos, second_pos;
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

bit_vector makeBP(vector<int>& array) {
    bit_vector bp;
    bp.resize(array.size()*2+2);
    int pos = 0;
    stack<int> Q;
    bp[pos] = 1;
    pos++;
    bp[pos] = 1;
    pos++;
    Q.push(array[0]);
    for(int i = 1;i < array.size();i++){
        while((not Q.empty()) and Q.top() < array[i]){
            bp[pos]=0;
            pos++;
            Q.pop();
        }
        Q.push(array[i]);
        bp[pos]=1;
        pos++;
    }
    while(not Q.empty()){
        bp[pos] = 0;
        pos++;
        Q.pop();
    }
    bp[pos]=0;
    pos++;
    return bp;
}
bit_vector makeSpine(vector<int> &array,bool reverse) {
    bit_vector spine;
    vector<vector<int>> leftSpineTable;
    vector<vector<int>> rightSpineTable;
    for(int i = 0;i<array.size();i++){
        vector<int> p;
        leftSpineTable.push_back(p);
        rightSpineTable.push_back(p);
    }
    stack<int> Q;
    for(int i = 0;i < array.size();i++) {
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
    int spine_size = 0;
    for(int i = 0;i < array.size();i++){
        int64_t leftsize = leftSpineTable[i].size();
        int64_t rightsize = rightSpineTable[i].size();
        spine_size += max(leftsize + rightsize - 1, (int64_t)0);
    }
    spine.resize(spine_size+2);
    int pos = 0;
    for(int i = 0;i < array.size();i++){
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
    return spine;
}


pair<int,int> getTop2FromSEA(vector<int>& A,rmq_succinct_rec<1024,128,true>& rmq,int startIndex,int endIndex){
    int r2m = -1;
    int rmqPos = rmq(startIndex,endIndex);
    if(rmqPos == startIndex){
        r2m = rmq(rmqPos+1,endIndex);
    }else if(rmqPos == endIndex){
        r2m = rmq(startIndex,rmqPos-1);
    }else{
        int left = rmq(startIndex,rmqPos-1);
        int right = rmq(rmqPos+1,endIndex);
        r2m = (A[left] < A[right] ? right : left);
    }
    return make_pair(rmqPos,r2m);
}
vector<pair<double,double>> dagTestUsingOptimization(vector<int>& example)//TEST DAG PART, return bpe,time pair for original, huffman, lzw
{
    queue<pair<int, int>> q;
    set<pair<int, int>> s;
    vector<bool> result;
    vector<pair<double, double>> testResult;
    bit_vector BP = makeBP(example);
    auto start_encoding = chrono::high_resolution_clock::now();
    rmq_succinct_rec<1024,128,true> rmq(&example);
    q.push(make_pair(0, example.size() - 1));
    s.insert(make_pair(0, example.size() - 1));
    while (!q.empty()) {
        bool isleaf;
        pair<int, int> range = q.front();
        q.pop();
        if(range.second - range.first < 2)
            continue;
        pair<int,int> top2_pos = getTop2FromSEA(example,rmq,range.first,range.second);
        int rmq_pos = top2_pos.first;
        bool rmq_side = (rmq_pos == range.first) || (rmq_pos == range.second) ? true : false;        
        int first_range_end = top2_pos.first > top2_pos.second ? top2_pos.first-1 : top2_pos.second-1;
        int second_range_start = top2_pos.first > top2_pos.second ? top2_pos.second+1 : top2_pos.first+1;
        pair<int, int> candidate_range_first = make_pair(range.first, first_range_end);
        pair<int, int> candidate_range_second = make_pair(second_range_start, range.second);
        if(s.find(candidate_range_first) == s.end()){
            s.insert(candidate_range_first);
            q.push(candidate_range_first);
        }
        if(s.find(candidate_range_second) == s.end()){
            s.insert(candidate_range_second);
            q.push(candidate_range_second);
        }
        if(!rmq_side){
            if (top2_pos.first > top2_pos.second) {//rmq on left
                result.push_back(false);
            } else {
                result.push_back(true);
            }
        }        
    }
    bit_vector result_bitvector = tobitvector(result);    
    auto end_encoding = chrono::high_resolution_clock::now();
    double timeEncoding = std::chrono::duration_cast<std::chrono::nanoseconds>(end_encoding - start_encoding).count() /(double) 1000;
    int huffDAGSize = getHuffmanSize(&result_bitvector);
    int huffBPSize = getHuffmanSize(&BP);    
    auto end_huffman = chrono::high_resolution_clock::now();
    double timeHuffman = std::chrono::duration_cast<std::chrono::nanoseconds>(end_huffman - end_encoding).count() /(double) 1000;
    int LZWDAGSize = getLz77Size(&result_bitvector);
    int LZWBPSize = getLz77Size(&BP);
    auto end_LZW = chrono::high_resolution_clock::now();
    double timeLZW = std::chrono::duration_cast<std::chrono::nanoseconds>(end_LZW - end_huffman).count() / (double)1000;
    double originalBPE = result_bitvector.size()/double(example.size()) + 2;
    double huffmanBPE = (huffDAGSize + huffBPSize)/double(example.size());
    double LZWBPE = (LZWDAGSize + LZWBPSize)/double(example.size());
    testResult.push_back(make_pair(originalBPE, timeEncoding));
    testResult.push_back(make_pair(huffmanBPE, timeHuffman));
    testResult.push_back(make_pair(LZWBPE, timeLZW));
    return testResult;
}
vector<pair<double,double>> dagTest(vector<int>& example)//TEST DAG PART, return bpe,time pair for original, huffman, lzw
{
    queue<pair<int, int>> q;
    set<pair<int, int>> s;
    vector<bool> result;
    vector<pair<double, double>> testResult;

    auto start_encoding = chrono::high_resolution_clock::now();
    bit_vector BP = makeBP(example);
    q.push(make_pair(0, example.size() - 1));
    s.insert(make_pair(0, example.size() - 1));
    while (!q.empty()) {
        bool isleaf;
        pair<int, int> range = q.front();
        q.pop();
        if(range.second - range.first < 2)
            continue;
        pair<int, int> top2_pos = second_largest_position_rmq_first(example, range.first, range.second);
        int rmq_pos = top2_pos.first;
        bool rmq_side = (rmq_pos == range.first) || (rmq_pos == range.second) ? true : false;        
        int first_range_end = top2_pos.first > top2_pos.second ? top2_pos.first-1 : top2_pos.second-1;
        int second_range_start = top2_pos.first > top2_pos.second ? top2_pos.second+1 : top2_pos.first+1;
        pair<int, int> candidate_range_first = make_pair(range.first, first_range_end);
        pair<int, int> candidate_range_second = make_pair(second_range_start, range.second);
        if(s.find(candidate_range_first) == s.end()){
            s.insert(candidate_range_first);
            q.push(candidate_range_first);
        }
        if(s.find(candidate_range_second) == s.end()){
            s.insert(candidate_range_second);
            q.push(candidate_range_second);
        }
        if(!rmq_side){
            if (top2_pos.first > top2_pos.second) {//rmq on left
                result.push_back(false);
            } else {
                result.push_back(true);
            }
        }        
    }
    bit_vector result_bitvector = tobitvector(result);    
    auto end_encoding = chrono::high_resolution_clock::now();
    double timeEncoding = std::chrono::duration_cast<std::chrono::nanoseconds>(end_encoding - start_encoding).count() /(double) 1000;
    int huffDAGSize = getHuffmanSize(&result_bitvector);
    int huffBPSize = getHuffmanSize(&BP);    
    auto end_huffman = chrono::high_resolution_clock::now();
    double timeHuffman = std::chrono::duration_cast<std::chrono::nanoseconds>(end_huffman - end_encoding).count() /(double) 1000;
    int LZWDAGSize = getLz77Size(&result_bitvector);
    int LZWBPSize = getLz77Size(&BP);
    auto end_LZW = chrono::high_resolution_clock::now();
    double timeLZW = std::chrono::duration_cast<std::chrono::nanoseconds>(end_LZW - end_huffman).count() / (double)1000;
    double originalBPE = result_bitvector.size()/double(example.size()) + 2;
    double huffmanBPE = (huffDAGSize + huffBPSize)/double(example.size());
    double LZWBPE = (LZWDAGSize + LZWBPSize)/double(example.size());
    testResult.push_back(make_pair(originalBPE, timeEncoding));
    testResult.push_back(make_pair(huffmanBPE, timeHuffman));
    testResult.push_back(make_pair(LZWBPE, timeLZW));
    return testResult;
}
vector<pair<double,double>> spineTest(vector<int>& example,bool reverse) {
    vector<pair<double, double>> testResult;

    auto start_encoding = chrono::high_resolution_clock::now();
    bit_vector BP = makeBP(example);
    bit_vector SPINE = makeSpine(example,reverse);
    auto end_encoding = chrono::high_resolution_clock::now();
    double timeEncoding = std::chrono::duration_cast<std::chrono::nanoseconds>(end_encoding - start_encoding).count() /(double) 1000;
    int huffBPSize = getHuffmanSize(&BP);
    int huffSPINESize = getHuffmanSize(&SPINE);
    auto end_huffman = chrono::high_resolution_clock::now();
    double timeHuffman = std::chrono::duration_cast<std::chrono::nanoseconds>(end_huffman - end_encoding).count() /(double) 1000;
    int LZWBPSize = getLz77Size(&BP);
    int LZWSPINESize = getLz77Size(&SPINE);
    auto end_LZW = chrono::high_resolution_clock::now();
    double timeLZW = std::chrono::duration_cast<std::chrono::nanoseconds>(end_LZW - end_huffman).count() / (double)1000;
    double originalBPE = (BP.size()+SPINE.size())/double(example.size());
    double huffmanBPE = (huffSPINESize + huffBPSize)/double(example.size());
    double LZWBPE = (LZWSPINESize + LZWBPSize)/double(example.size());
    testResult.push_back(make_pair(originalBPE, timeEncoding));
    testResult.push_back(make_pair(huffmanBPE, timeHuffman));
    testResult.push_back(make_pair(LZWBPE, timeLZW));

    return testResult;
}