#include <iostream>
#include <fstream> 
#include <vector>
#include <string>
#include <algorithm>
#include <limits>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <functional>

#define idx(i, j, w) ((i) * (w) + (j))



__global__
void gotohAUX2(int diag, int m, int n, const char* A,
    const char* B,
    int openGap,
    int extendGap, const int* submat){
          size_t index = blockIdx.x * blockDim.x + threadIdx.x;
          int i = max(1, diag-n) + index;
          int j = diag - i;


    }

void gotoch_align_cuda(const string&A, const string &B, int openGap, int extendHap,  const vector<vector<int>> &submat) {
    
    const int THREADS_PER_BLOCK = 128;
    
    
    const int NEG_INF = numeric_limits<int>::min()/2;
    const int m = A.size();
    const int n = B.size();
    const int width = n + 1;
    const int size = (m + 1) * (n + 1);


    //here we are allocating host flattened matrices
    vector<int> hM(size, NEG_INF);
    vector<int> hI(size, NEG_INF);
    vector<int> hD(size, NEG_INF);

}
  

    

    //mobving data to device

