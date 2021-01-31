#include <iostream>
#include <vector>
#include <stdint.h>
#include "DFUDSEncoding.hpp"
using namespace std;
int main()
{
    vector<uint> A = {2,10,3,0,11,1,8,6,7,9,4,5};
    DFUDSEncoding DFUDSencoding(A);
    cout<<"DFUDS SEQUENCE"<<endl;
    for(int i=0;i<DFUDSencoding.dfuds.size();i++)
        cout<<DFUDSencoding.dfuds[i];
    cout<<endl;
    cout<<"SPINE SEQUENCE"<<endl;
    for(int i=0;i<DFUDSencoding.spine.size();i++)
        cout<<DFUDSencoding.spine[i];
    cout<<endl;
    cout<<"RT2Q QUERY RESULT R2TQ(3,9)"<<endl;
    auto result = DFUDSencoding.top2_pos(2,8);
    cout<<"RMQ INDEX : "<<result.first<<" , TOP2 : "<<result.second<<endl;

}