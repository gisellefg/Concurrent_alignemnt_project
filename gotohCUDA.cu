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
#include <cuda_runtime.h>
#include <cuda.h>
using namespace std;

#define IDX(i, j, w) ((i) * (w) + (j))


//this function is a function to read fasta file, i used GPT to write it
string fastaReader(const string &path) {
    ifstream in(path);
    if (!in) throw runtime_error("Cannot open FASTA file: " + path);
    string line, seq;
    while (getline(in, line)) {
        if (line.empty() || line[0] == '>') continue;
        seq += line;
    }
    return seq;
}

//For comparing results 
struct AlignmentResult {
    std::string A_aligned;
    std::string B_aligned;
    int score;

    AlignmentResult(const std::string& a, const std::string& b, int s)
        : A_aligned(a), B_aligned(b), score(s) {}


};

// for transferring variables between function calls
struct AlignmentMatrices {
    std::vector<int> M, I, D;
    std::vector<int> traceM, traceI, traceD;
    int m, n, width;
};



//Computes one anti-diagonal
__global__ void gotohAUX2(int diag, int m, int n, const char* A, const char* B,
                int openGap, int extendGap, const int* submat,
                int* M, int* I, int* D, int*traceM, int* traceI, int*traceD){
            size_t index = blockIdx.x * blockDim.x + threadIdx.x;
            int i = max(1, diag-n) + index;
            int j = diag - i;

            if (i > m || j < 1 || j > n) return;

            int width = n + 1;
            int idx       = IDX(i, j, width);
            int up        = IDX(i - 1, j, width);
            int left      = IDX(i, j - 1, width);
            int diagPrev  = IDX(i - 1, j - 1, width);



            int sub = submat[(int)A[i - 1] * 128 + (int)B[j - 1]];

            // Compute I (insertion score)
            //Here we open new gap vertically
            int openI = M[up] - (openGap + extendGap);
            //Extend existing vertical gap
            int extI  = I[up] - extendGap;
            if (openI >= extI) {
                I[idx] = openI;
                traceI[idx] = 0; 
            } else {
                I[idx] = extI;
                traceI[idx] = 1;
            }

            
            // Compute D (deletion score)
            int openD = M[left] - (openGap + extendGap);
            int extD  = D[left] - extendGap;
            if (openD >= extD) {
                D[idx] = openD;
                traceD[idx] = 0;
            } else {
                D[idx] = extD;
                traceD[idx] = 1;
            }



            // Compute M (match/mismatch score)
            //We have three possibilities here: match/mismatch from M, gap from I or gap from D
            int diagS = M[diagPrev] + sub;

            if (diagS >= I[idx] && diagS >= D[idx]) {
                M[idx] = M[diagPrev] + sub;
                traceM[idx] = 0;

            } else if (I[idx] >= D[idx]) {
                M[idx] = I[idx];
                traceM[idx] = 1;

            } else {
                M[idx] = D[idx];
                traceM[idx] = 2;
            }

    }

AlignmentMatrices gotoch_align_cuda(const string&A, const string &B, int openGap, int extendGap, const std::vector<std::vector<int>> &submat) {
    
    const int THREADS_PER_BLOCK = 128;
    const int NEG_INF = numeric_limits<int>::min()/2;
    
    char *d_A, *d_B;
    int *d_M, *d_I, *d_D, *d_submat;
    int *d_traceM, *d_traceI, *d_traceD;

    

    
    const int m = A.size();
    const int n = B.size();
    const int width = n + 1;
    const int size = (m + 1) * (n + 1);

    //Flatten substition matrix (code received from CHATGPT)
    std::vector<int> flatSubmat(128 * 128, 0);
    for (int i = 0; i < 128; ++i)
        for (int j = 0; j < 128; ++j)
            flatSubmat[i * 128 + j] = submat[i][j];

    
    

    //here we are allocating host flattened matrices (HOST)
    vector<int> hM(size, NEG_INF);
    vector<int> hI(size, NEG_INF);
    vector<int> hD(size, NEG_INF);
    hM[0] = hI[0] = hD[0] = 0;

    for (int i = 1; i <= m; i++) {
        hI[i * width] =  -(openGap +(i-1)*extendGap);
    }
    for (int j = 1; j <= n; j++) {
        hD[j] = -(openGap + (j-1)*extendGap);
    }


    vector<int> hTraceM(size);
    vector<int> hTraceI(size);
    vector<int> hTraceD(size);




    // Allocate to host
    cudaMalloc(&d_A, m*sizeof(char));
    cudaMalloc(&d_B, n*sizeof(char));
    cudaMemcpy(d_A, A.c_str(), m, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, B.c_str(), n, cudaMemcpyHostToDevice);

    cudaMalloc(&d_M, size * sizeof(int));
    cudaMalloc(&d_I, size * sizeof(int));
    cudaMalloc(&d_D, size * sizeof(int));
    cudaMemcpy(d_M, hM.data(), size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_I, hI.data(), size * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_D, hD.data(), size * sizeof(int), cudaMemcpyHostToDevice);

    cudaMalloc(&d_submat, 128 * 128 * sizeof(int));
    cudaMemcpy(d_submat, flatSubmat.data(), 128 * 128 * sizeof(int), cudaMemcpyHostToDevice); 


    //traceback matrices
    cudaMalloc(&d_traceM, size * sizeof(int));
    cudaMalloc(&d_traceI, size * sizeof(int));
    cudaMalloc(&d_traceD, size * sizeof(int));

    //Compute on GPU, we launch kernels for each diagonal
    for (int diag = 2; diag <= m + n; ++diag){
        int i_min = max(1, diag - n);
        int i_max = min(m, diag - 1);
        int total = i_max - i_min + 1;
        if (total <= 0) continue;

        int blocks = (total + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    
        gotohAUX2<<<blocks, THREADS_PER_BLOCK>>>(diag, m, n, d_A, d_B,
                                         openGap, extendGap, d_submat,
                                         d_M, d_I, d_D,
                                         d_traceM, d_traceI, d_traceD);


        cudaDeviceSynchronize();


    }

    //Results back to host
    cudaMemcpy(hM.data(), d_M, size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(hI.data(), d_I, size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(hD.data(), d_D, size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(hTraceM.data(), d_traceM, size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(hTraceI.data(), d_traceI, size * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(hTraceD.data(), d_traceD, size * sizeof(int), cudaMemcpyDeviceToHost);


    // compare values of different matrices




    //Free CUDA memory
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_M);
    cudaFree(d_I);
    cudaFree(d_D);
    cudaFree(d_traceM);
    cudaFree(d_traceI);
    cudaFree(d_traceD);
    cudaFree(d_submat);

    return {
        std::move(hM), std::move(hI), std::move(hD),
        std::move(hTraceM), std::move(hTraceI), std::move(hTraceD),
        m, n, width
    };

}

AlignmentResult gotoch_align_result(const std::string& A, const std::string& B,
                                const std::vector<int>& hM,
                                const std::vector<int>& hI,
                                const std::vector<int>& hD,
                                const std::vector<int>& hTraceM,
                                const std::vector<int>& hTraceI,
                                const std::vector<int>& hTraceD,
                                int m, int n, int width) {
    auto IDX = [width](int i, int j) { return i * width + j; };


    int end = IDX(m, n, width);
    int score_D = hD[end];
    int score_M = hM[end];
    int score_I = hI[end];
    int score;
    char chosen;

    if (score_M >= score_I && score_M >= score_D) {
        score = score_M;
        chosen = 'M';
    } else if (score_D >= score_I) {
        score = score_D;
        chosen = 'D';
    }else {
        score = score_I;
        chosen = 'I';
    }
    
    char current = chosen;

    //TO DO: add traceback logic as in CPU
    std::string A_aligned, B_aligned;

    int i = m, j = n;

    while (i > 0 || j > 0) {
        if (current == 'M') {
            if (i == 0) {
                current = 'D';
                continue;
            }
            if (j == 0) {
                current = 'I';
                continue;
            }
            int trace = hTraceM[IDX(i, j, width)];
            if (trace == 0) { // from M[i-1][j-1]
                A_aligned += A[i - 1];
                B_aligned += B[j - 1];
                i--; j--;
            } else if (trace == 1) {
                current = 'I';
            } else {
                current = 'D';
            }
        } else if (current == 'I') {
            if (i == 0) {
                current = 'D';
                continue;
            }
            A_aligned += A[i - 1];
            B_aligned += '-';
            if (hTraceI[IDX(i, j, width)] == 0) {
                current = 'M';
            } else {
                current = 'I';
            }
            i--;
        } else { // current == 'D'
            if (j == 0) {
                current = 'I';
                continue;
            }
            B_aligned += B[j - 1];
            A_aligned += '-';
            if (hTraceD[IDX(i, j, width)] == 0) {
                current = 'M';
            } else {
                current = 'D';
            }
            j--;
        }
    }

    //we reverse so that the strings we get are from left to right 
    std::reverse(A_aligned.begin(), A_aligned.end());
    std::reverse(B_aligned.begin(), B_aligned.end());

    return AlignmentResult(A_aligned, B_aligned, score);
}


int main() {

    //testing code with input strings
    std::string A = fastaReader("seqA.fasta");
    std::string B = fastaReader("seqB.fasta");

    // parameters for scoring alignment
    int openGap = 10;
    int extendGap = 1;
    int match = 3;
    int mismatch = -1;

    // submatrix initialization
    vector<vector<int>> submat(128, vector<int>(128, mismatch));
    for (char c : {'A','C','G','T'}) submat[c][c] = match;

    // get results (body of code)
    AlignmentMatrices ms = gotoch_align_cuda(A, B, openGap, extendGap, submat);
    AlignmentResult result = gotoch_align_result(A, B, ms.M, ms.I, ms.D, ms.traceM, ms.traceI, ms.traceD, ms.m, ms.n, ms.width);

    std::cout << "Aligned A: " << result.A_aligned << "\n";
    std::cout << "Aligned B: " << result.B_aligned << "\n";
    std::cout << "Alignment score: " << result.score << "\n";

    return 0;
}
  

    


  

    

