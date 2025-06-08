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
#include "fasta_reader.h"
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
 

std::tuple<std::vector<int>, std::vector<int>, std::vector<int>> compute_linear_space_cpu(
    const std::string& A, const std::string& B, int openGap, int extendGap,
    const std::vector<std::vector<int>>& submat, bool reverse = false) {

    int m = A.size();
    int n = B.size();
    int w = n + 1;
    const int NEG_INF = std::numeric_limits<int>::min() / 2;

    //Only two rows will be kept in memory at any time to optimize memory complexity
    std::vector<int> hM0(w), hM1(w);
    std::vector<int> hI0(w), hI1(w);
    std::vector<int> hD0(w), hD1(w);

    hM0[0] = 0;
    hI0[0] = NEG_INF;
    hD0[0] = NEG_INF;

    for (int j = 1; j <=n; ++j) {
        hM0[j] = NEG_INF;
        hD0[j] = -openGap - (j - 1) * extendGap;
        hI0[j] = NEG_INF;
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

        //Swap current and previous rows
        std::swap(hM0, hM1);
        std::swap(hI0, hI1);
        std::swap(hD0, hD1);
    }

    return {hM0, hI0, hD0};
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


/*
int main() {
    std::string fastaFile = "uniref50.fasta";
    int length = 5000;
    int openGap = 5;
    int extendGap = 1;
    int match = 3;
    int mismatch = -1;

    // Compare sequences from [x, y) vs [x1, y1)
    int x = 0, y = 5;     // first range
    int x1 = 5, y1 = 10;  // second range

    auto sequences = loadAllSequences(fastaFile);
    if (sequences.empty()) {
        std::cerr << "No sequences loaded.\n";
        return 1;
    }

    std::vector<std::vector<int>> submat(128, std::vector<int>(128, mismatch));
    for (char c : {'A', 'C', 'G', 'T'}) submat[c][c] = match;

    for (int i = x; i < y && i < sequences.size(); ++i) {
        std::string A = sequences[i].substr(0, length);

        for (int j = x1; j < y1 && j < sequences.size(); ++j) {
            std::string B = sequences[j].substr(0, length);
            auto [finalM, finalI, finalD] = compute_linear_space_cpu(A, B, openGap, extendGap, submat);
            int score = std::max({finalM.back(), finalI.back(), finalD.back()});
            std::cout << i << " vs " << j << " → score: " << score << "\n";
        }
    }

    return 0;
}


*/


int main() {
    std::string fastaFile = "uniref50.fasta";
    int openGap = 5;
    int extendGap = 1;
    int match = 3;
    int mismatch = -1;

    // Sequence index ranges (inclusive)
    int startA = 1, endA = 5;
    int startB = 6, endB = 10;

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
    //Here we compare two concatenated sequences
    std::string A, B;
    for (int i = startA; i <= endA; ++i) A += sequences[i];
    for (int i = startB; i <= endB; ++i) B += sequences[i];

    std::vector<std::vector<int>> submat(128, std::vector<int>(128, mismatch));
    for (char c : {'A', 'C', 'G', 'T'}) submat[c][c] = match;

    auto [finalM, finalI, finalD] = compute_linear_space_cpu(A, B, openGap, extendGap, submat);
    int score = std::max({finalM.back(), finalI.back(), finalD.back()});
    std::cout << "Alignment score (seq " << startA << "-" << endA
              << " vs " << startB << "-" << endB << "): " << score << "\n";

    return 0;
}


    

