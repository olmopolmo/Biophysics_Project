#include<iostream>
#include<vector>
#include<cmath>
#include "graph.hpp"
using namespace std;

int main(){
    int n;
    int i;
    int j;
    int distance;
    int x;
    cin >> n >> i >> j;
    for (int a=0; a<n ; a++){
        distance = 0;
        for (int b=0; b<j;b++){
            cin >> x;
            distance += x;
        }
    cout << distance << endl;
    }
}