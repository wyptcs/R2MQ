#include <iostream>
#include <vector>
#include <algorithm>
#include <bits/stdc++.h>
#include <stack>
#include <queue>
#include <string>
#include <unordered_map>
#include <array>
#include <chrono>
#include "huffman.cpp"
#include <sdsl/rrr_vector.hpp>
using namespace std;
using namespace sdsl;
pair<int,int> second_largest_position_rmq_first(vector<int>& arr,int start,int end);
bit_vector tobitvector(vector<bool>& input);
pair<double,int> dag_test(vector<int>& example)//TEST DAG PART, RETURN CONSTRUCTION TIME, EXCLUDE HUFFMAN
{
    queue<pair<int, int>> q;
    set<pair<int, int>> s;
    vector<bool> result;
    //pair<int,int> second_largest_pos = second_largest_position(example,0,example.size()-1);
    auto start_encoding = chrono::high_resolution_clock::now();
    q.push(make_pair(0, example.size() - 1));
    s.insert(make_pair(0, example.size() - 1));
    while (!q.empty()) {
        bool isleaf;
        pair<int, int> range = q.front();
        pair<int, int> top2_pos = second_largest_position_rmq_first(example, range.first, range.second);
        int rmq_pos = example[top2_pos.first] < example[top2_pos.second] ? top2_pos.first : top2_pos.second;
        bool rmq_side = (rmq_pos == range.first) || (rmq_pos == range.second) ? true : false;

        if (range.second - range.first < 3)
            isleaf = true;
        else
            isleaf = false;
        q.pop();
        int first_range_end = top2_pos.second - 1;
        int second_range_start = top2_pos.first + 1;
        pair<int, int> candidate_range_first = make_pair(range.first, first_range_end);
        pair<int, int> candidate_range_second = make_pair(second_range_start, range.second);

        if (!isleaf && s.find(candidate_range_first) == s.end() && first_range_end - range.first > 1) {
            pair<int,int> top2_pos_candidate = second_largest_position_rmq_first(example,candidate_range_first.first,candidate_range_first.second);
            rmq_pos = example[top2_pos_candidate.first] < example[top2_pos_candidate.second] ? top2_pos_candidate.first : top2_pos_candidate.second;
            rmq_side = (rmq_pos == candidate_range_first.first) || (rmq_pos == candidate_range_first.second) ? true : false;

            if(!rmq_side) {
                s.insert(candidate_range_first);
                q.push(candidate_range_first);
            }
        }
        if (!isleaf && s.find(candidate_range_second) == s.end() && range.second - second_range_start > 1) {
            pair<int,int> top2_pos_candidate = second_largest_position_rmq_first(example,candidate_range_second.first,candidate_range_second.second);
            rmq_pos = example[top2_pos_candidate.first] < example[top2_pos_candidate.second] ? top2_pos_candidate.first : top2_pos_candidate.second;
            rmq_side = (rmq_pos == candidate_range_second.first) || (rmq_pos == candidate_range_second.second) ? true : false;

            if(!rmq_side) {
                s.insert(candidate_range_second);
                q.push(candidate_range_second);
            }
        }
        if (example[top2_pos.first] > example[top2_pos.second]) {//rmq on left
            result.push_back(false);
        } else {
            result.push_back(true);
        }
    }
    bit_vector result_bitvector = tobitvector(result);
    auto end_encoding = chrono::high_resolution_clock::now();
    double time = std::chrono::duration_cast<std::chrono::nanoseconds>(end_encoding - start_encoding).count() /(double) 1000;
    int huffSize = getHuffmanSize(&result_bitvector);
    return make_pair(time,huffSize);
    //return make_pair(result.size(),huffSize);
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
bit_vector tobitvector(vector<bool>& input){
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