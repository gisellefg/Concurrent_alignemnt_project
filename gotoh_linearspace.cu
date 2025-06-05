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
#include "gotohCUDA.h"
#include "types.h"

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
    
__global__ void gotohAUX2(int i, int n, const char*A, const char*B, int openGap, 
int extendGap, const int*submat, const int*pM, const int* pI, const int* pD, int* cM, int*cI, int* cD){
            
            //index at 0 is handled on CPU 
            size_t index = blockIdx.x * blockDim.x + threadIdx.x + 1;
            int sub = submat[(int)A[i - 1] * 128 + (int)B[index - 1]];
            if (index > n){
                return;
            }

            // Compute I (insertion score)
            //Here we open new gap vertically
            int openI = pM[index] - (openGap + extendGap);
            //Extend existing vertical gap
            int extI  = pI[index] - extendGap;
            cI[index] = max(openI, extI);

            
            // Compute D (deletion score)
            int openD = cM[index - 1] - (openGap + extendGap);
            int extD  = cD[index - 1] - extendGap;
            cD[index] = max(openD, extD);


            // Compute M (match/mismatch score)
            int diag = pM[index-1] + sub;
            cM[index] = max(max(cI[index], cD[index]), diag);


};

std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> compute_linear_space_cpu(
    const std::string& A, const std::string& B, int openGap, int extendGap,
    const std::vector<std::vector<int>>& submat, bool reverse = false) {

    int m = A.size();
    int n = B.size();
    int w = n + 1;
    const int NEG_INF = std::numeric_limits<int>::min() / 2;

    std::vector<int> hM0(w, NEG_INF), hM1(w, NEG_INF);
    std::vector<int> hI0(w, NEG_INF), hI1(w, NEG_INF);
    std::vector<int> hD0(w, NEG_INF), hD1(w, NEG_INF);

    hM0[0] = hI0[0] = hD0[0] = 0;
    for (int j = 1; j < w; ++j) {
        hM0[j] = NEG_INF;
        hI0[j] = -openGap - (j - 1) * extendGap;
        hD0[j] = NEG_INF;
    }

    for (int i = 1; i <= m; ++i) {
        hM1[0] = hD1[0] = NEG_INF;
        hI1[0] = -openGap - (i - 1) * extendGap;

        for (int j = 1; j <= n; ++j) {
            int sub = submat[(int)A[i - 1]][(int)B[j - 1]];

            // Insertions
            int openI = hM0[j] - (openGap + extendGap);
            int extI  = hI0[j] - extendGap;
            hI1[j] = std::max(openI, extI);

            // Deletions
            int openD = hM1[j - 1] - (openGap + extendGap);
            int extD  = hD1[j - 1] - extendGap;
            hD1[j] = std::max(openD, extD);

            // Match/Mismatch
            int diag = hM0[j - 1] + sub;
            hM1[j] = std::max({hI1[j], hD1[j], diag});
        }

        std::swap(hM0, hM1);
        std::swap(hI0, hI1);
        std::swap(hD0, hD1);
    }

    return {hM0, hI0, hD0};
}

/*
//we want to return all three final rows
tuple<vector<int>, vector<int>, vector<int>> compute_linear_space(
    const string& A, const string& B, int openGap, int extendGap, 
    const vector<vector<int>>& submat, bool reverse = false) {

    int m = A.size();
    int n = B.size();
    int w = n + 1;
    const int THREADS_PER_BLOCK = 128;
    const int NEG_INF = numeric_limits<int>::min() / 2;

    vector<int> hM0(w, NEG_INF), hM1(w, NEG_INF);
    vector<int> hI0(w, NEG_INF), hI1(w, NEG_INF);
    vector<int> hD0(w, NEG_INF), hD1(w, NEG_INF);
    hM0[0] = hI0[0] = hD0[0] = 0;
    for (int j = 1; j < w; ++j) {
        hM0[j] = NEG_INF;
        hI0[j] = -openGap - (j - 1) * extendGap;
        hD0[j] = NEG_INF;
    }

    vector<int> flatSubmat(128 * 128, 0);
    for (int i = 0; i < 128; ++i)
        for (int j = 0; j < 128; ++j)
            flatSubmat[i * 128 + j] = submat[i][j];

    char *dA, *dB;
    int *dM0, *dM1, *dI0, *dI1, *dD0, *dD1, *dSubmat;

    cudaMalloc(&dA, m);
    cudaMalloc(&dB, n);
    cudaMemcpy(dA, A.data(), m, cudaMemcpyHostToDevice);
    cudaMemcpy(dB, B.data(), n, cudaMemcpyHostToDevice);

    cudaMalloc(&dM0, w * sizeof(int)); cudaMalloc(&dM1, w * sizeof(int));
    cudaMalloc(&dI0, w * sizeof(int)); cudaMalloc(&dI1, w * sizeof(int));
    cudaMalloc(&dD0, w * sizeof(int)); cudaMalloc(&dD1, w * sizeof(int));

    cudaMemcpy(dM0, hM0.data(), w * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dI0, hI0.data(), w * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dD0, hD0.data(), w * sizeof(int), cudaMemcpyHostToDevice);

    cudaMalloc(&dSubmat, 128 * 128 * sizeof(int));
    cudaMemcpy(dSubmat, flatSubmat.data(), 128 * 128 * sizeof(int), cudaMemcpyHostToDevice);

    int *pM = dM0, *pI = dI0, *pD = dD0;
    int *cM = dM1, *cI = dI1, *cD = dD1;

    for (int i = 1; i <= m; ++i) {
        hI1[0] = -openGap - i * extendGap;
        hM1[0] = hD1[0] = NEG_INF;

        int blocks = (n + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
        gotohAUX2<<<blocks, THREADS_PER_BLOCK>>>(i, n, dA, dB,
                                                 openGap, extendGap, dSubmat,
                                                 pM, pI, pD,
                                                 cM, cI, cD);

        std::swap(pM, cM);
        std::swap(pI, cI);
        std::swap(pD, cD);
    }

    cudaDeviceSynchronize();

    vector<int> final_rowM(w), final_rowI(w), final_rowD(w);
    cudaMemcpy(final_rowM.data(), pM, w * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(final_rowI.data(), pI, w * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(final_rowD.data(), pD, w * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(dA); cudaFree(dB); cudaFree(dSubmat);
    cudaFree(dM0); cudaFree(dM1);
    cudaFree(dI0); cudaFree(dI1);
    cudaFree(dD0); cudaFree(dD1);

    return {final_rowM, final_rowI, final_rowD};
}

*/





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


ScoreTime alignCPU(const std::string& A, const std::string& B,
            const int openGap, const int extendGap,
            const int match, const int mismatch) {

    vector<vector<int>> submat(128, vector<int>(128, mismatch)); // initialize submatrix
    for (char c : {'A','C','G','T'}) submat[c][c] = match;

    cudaEvent_t cpu_start, cpu_end;
    cudaEventCreate(&cpu_start);
    cudaEventCreate(&cpu_end);

    cudaEventRecord(cpu_start); // start time
    // In aligncpu()
    auto [finalM, finalI, finalD] = compute_linear_space_cpu(A, B, openGap, extendGap, submat); 
    cudaEventRecord(cpu_end); // end time
    cudaEventSynchronize(cpu_end);

    float cpu_time_ms;
    cudaEventElapsedTime(&cpu_time_ms, cpu_start, cpu_end);

    cudaEventDestroy(cpu_start);
    cudaEventDestroy(cpu_end);


    int score = std::max({finalM.back(), finalI.back(), finalD.back()});
    return ScoreTime(score, cpu_time_ms);
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
    std::vector<vector<int>> submat(128, vector<int>(128, mismatch));
    for (char c : {'A','C','G','T'}) submat[c][c] = match;

    // timing
    cudaEvent_t cpu_start, cpu_end;
    cudaEventCreate(&cpu_start);
    cudaEventCreate(&cpu_end);

    //Since we no longer need to calculate the traceback, we only need the final row vector
    cudaEventRecord(cpu_start);
    auto [finalM, finalI, finalD] = compute_linear_space_cpu(A, B, openGap, extendGap, submat);
    cudaEventRecord(cpu_end);
    cudaEventSynchronize(cpu_end);

    float cpu_time_ms;
    cudaEventElapsedTime(&cpu_time_ms, cpu_start, cpu_end);
    std::cout << "cpu Time: " << cpu_time_ms / 1000.0 << " s\n";

    std::cout << "M[m][n] = " << finalM.back() << "\n";
    std::cout << "I[m][n] = " << finalI.back() << "\n";
    std::cout << "D[m][n] = " << finalD.back() << "\n";

    // get results
    int score = std::max({finalM.back(), finalI.back(), finalD.back()});
    std::cout << "Alignment score: " << score << "\n";

    cudaEventDestroy(cpu_start);
    cudaEventDestroy(cpu_end);


    return 0;
}

    



  

    

