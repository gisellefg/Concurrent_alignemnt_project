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
#include <tuple>
#include "gotohCUDA.h"
#include "types.h"
using namespace std;
#include <cuda_runtime.h>
#include <cuda.h>
#include "fasta_reader.h"

#define IDX(i, j, w) ((i) * (w) + (j))


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
 

__global__ void gotohAUX2(
    int i, int n, const char* A, const char* B, const int* submat,
    int openGap, int extendGap,
    const int* pM, const int* pI, const int* pD,
    int* cM, int* cI, int* cD
) {
    int j = blockIdx.x * blockDim.x + threadIdx.x + 1;
    if (j > n) return;

    int sub = submat[(int)A[i - 1] * 128 + (int)B[j - 1]];

    int openI = pM[j] - (openGap + extendGap);
    int extI  = pI[j] - extendGap;
    cI[j] = max(openI, extI);

    int openD = cM[j - 1] - (openGap + extendGap);
    int extD  = cD[j - 1] - extendGap;
    cD[j] = max(openD, extD);

    int diag = pM[j - 1] + sub;
    cM[j] = max(max(cI[j], cD[j]), diag);
}


__host__ tuple<vector<int>, vector<int>, vector<int>> compute_linear_space_gpu(
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

    cudaMalloc(&dM0, w * sizeof(int)); 
    cudaMalloc(&dM1, w * sizeof(int));
    cudaMalloc(&dI0, w * sizeof(int)); 
    cudaMalloc(&dI1, w * sizeof(int));
    cudaMalloc(&dD0, w * sizeof(int)); 
    cudaMalloc(&dD1, w * sizeof(int));

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
        gotohAUX2<<<blocks, THREADS_PER_BLOCK>>>(i, n, dA, dB, dSubmat,
        openGap, extendGap,
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

    cudaFree(dA); 
    cudaFree(dB); 
    cudaFree(dSubmat);
    cudaFree(dM0); 
    cudaFree(dM1);
    cudaFree(dI0); 
    cudaFree(dI1);
    cudaFree(dD0); 
    cudaFree(dD1);

    return {final_rowM, final_rowI, final_rowD};
}

std::vector<std::string> loadAllSequences(const std::string& path) {
    std::ifstream in(path);
    std::vector<std::string> seqs;
    std::string line, seq;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!seq.empty()) { seqs.push_back(seq); seq.clear(); }
        } else {
            seq += line;
        }
    }
    if (!seq.empty()) seqs.push_back(seq);
    return seqs;
}


int main() {
    std::string fastaFile = "uniref50.fasta";
    int openGap = 5;
    int extendGap = 1;
    int match = 3;
    int mismatch = -1;

    // Sequence index ranges (inclusive)
    int startA = 1, endA = 5;
    int startB = 101, endB = 105;

    auto sequences = loadAllSequences(fastaFile);
    if (sequences.empty()) {
        std::cerr << "No sequences loaded.\n";
        return 1;
    }

    if (endA >= sequences.size() || endB >= sequences.size()) {
        std::cerr << "Index out of bounds. Total sequences loaded: " << sequences.size() << "\n";
        return 1;
    }

    // Concatenate sequences from startA to endA (inclusive)
    std::string A, B;
    for (int i = startA; i <= endA; ++i) A += sequences[i];
    for (int i = startB; i <= endB; ++i) B += sequences[i];

    std::vector<std::vector<int>> submat(128, std::vector<int>(128, mismatch));
    for (char c : {'A', 'C', 'G', 'T'}) submat[c][c] = match;

    auto [finalM, finalI, finalD] = compute_linear_space_gpu(A, B, openGap, extendGap, submat);
    int score = std::max({finalM.back(), finalI.back(), finalD.back()});
    std::cout << "Alignment score (seq " << startA << "-" << endA
              << " vs " << startB << "-" << endB << "): " << score << "\n";

    return 0;
}
