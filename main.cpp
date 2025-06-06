#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include "gotohCPU.h"
#include "gotohCUDA.h"
#include "gotoh_linearspace.h"
#include "types.h"
#include<chrono>
#include "fasta_reader.cpp"
#include <cuda_runtime.h>
#include "gotoh_linearspace.h"
using namespace std;

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


struct Strings {
    std::string A;
    std::string B;
};





// Batch Test functions are for computing computation times of large batches of sequences at once. Separate implementations for CPU and GPU.

void BatchTestCPU(const std::vector<Strings>& seqs, const int open, const int extend, const int match, const int mismatch) {
    for (size_t i = 0; i < seqs.size(); ++i) {
        ScoreTime scoretime = alignCPU(seqs[i].A, seqs[i].B, open, extend, match, mismatch);

        std::cout << "test:  " << i+1 << " score (CPU): " << scoretime.score 
                  << " time taken: " << scoretime.time << " ms" << std::endl;
    }
}

void BatchTestGPU(const std::vector<Strings>& seqs, const int open, const int extend, const int match, const int mismatch) {
    for (size_t i = 0; i < seqs.size(); ++i) {
        ScoreTime scoretime = alignGPU(seqs[i].A, seqs[i].B, open, extend, match, mismatch);

        std::cout << "test: " << i+1 << " score (GPU): " << scoretime.score 
                  << " time taken: " << scoretime.time << " ms" << std::endl;
    }
}

ScoreTime alignNW_Affine(const std::string& A, const std::string& B,
                                  int match, int mismatch, int openGap, int extendGap) {
    const int NEG_INF = numeric_limits<int>::min()/2;
    int m = A.size();
    int n = B.size();

    std::vector<std::vector<int>> M(m + 1, std::vector<int>(n + 1, NEG_INF));
    std::vector<std::vector<int>> Ix(m + 1, std::vector<int>(n + 1, NEG_INF)); // Gap in B (vertical)
    std::vector<std::vector<int>> Iy(m + 1, std::vector<int>(n + 1, NEG_INF)); // Gap in A (horizontal)

    // Initialization
    M[0][0] = 0;

    // Initialize first row and column
    // Cost for a gap of length k = openGap + k * extendGap
    for (int i = 1; i <= m; ++i) {
        // Align A[0...i-1] with i gaps in B
        int gap_penalty_i = openGap + i * extendGap;
        Ix[i][0] = -gap_penalty_i;
        M[i][0] = -gap_penalty_i; 
        Iy[i][0] = NEG_INF;       
    }
    for (int j = 1; j <= n; ++j) {
        // Align B[0...j-1] with j gaps in A
        int gap_penalty_j = openGap + j * extendGap;
        Iy[0][j] = -gap_penalty_j;
        M[0][j] = -gap_penalty_j; 
        Ix[0][j] = NEG_INF;       
    }
    
    auto start = std::chrono::high_resolution_clock::now();

    // Fill DP matrices
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int s = (A[i - 1] == B[j - 1]) ? match : mismatch;

            // Option 1: Open a new gap in B (transition from M[i-1][j])
            // Penalty for first char in gap: openGap + extendGap
            int open_ix = (M[i - 1][j] == NEG_INF) ? NEG_INF : M[i - 1][j] - (openGap + extendGap);
            // Option 2: Extend an existing gap in B (transition from Ix[i-1][j])
            // Penalty for extending: extendGap
            int extend_ix = (Ix[i - 1][j] == NEG_INF) ? NEG_INF : Ix[i - 1][j] - extendGap;
            Ix[i][j] = std::max(open_ix, extend_ix);

            // Option 1: Open a new gap in A (transition from M[i][j-1])
            // Penalty for first char in gap: openGap + extendGap
            int open_iy = (M[i][j - 1] == NEG_INF) ? NEG_INF : M[i][j - 1] - (openGap + extendGap);
            // Option 2: Extend an existing gap in A (transition from Iy[i][j-1])
            // Penalty for extending: extendGap
            int extend_iy = (Iy[i][j - 1] == NEG_INF) ? NEG_INF : Iy[i][j - 1] - extendGap;
            Iy[i][j] = std::max(open_iy, extend_iy);
            
            // M[i][j]: A[i-1] is aligned with B[j-1] (match/mismatch)
            int match_mismatch_score = (M[i - 1][j - 1] == NEG_INF) ? NEG_INF : M[i - 1][j - 1] + s;
            
            M[i][j] = std::max({match_mismatch_score, Ix[i][j], Iy[i][j]});
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    float time_ms = std::chrono::duration<float, std::milli>(end - start).count();

    return ScoreTime(M[m][n], time_ms);
}

ScoreTime alignlinearCPU(const std::string& A, const std::string& B,
                         int openGap, int extendGap, int match, int mismatch) {
    std::vector<std::vector<int>> submat(128, std::vector<int>(128, mismatch));
    for (char c : {'A', 'C', 'G', 'T'}) submat[c][c] = match;

    cudaEvent_t start, end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start);

    auto [M, I, D] = compute_linear_space_cpu(A, B, openGap, extendGap, submat);

    cudaEventRecord(end);
    cudaEventSynchronize(end);

    float time_ms;
    cudaEventElapsedTime(&time_ms, start, end);
    cudaEventDestroy(start);
    cudaEventDestroy(end);

    int score = std::max({M.back(), I.back(), D.back()});
    return ScoreTime(score, time_ms);
}


// compare the CPU and GPU run times for ONE pair of sequences (not a batch test)
void testCPUandGPU(const Strings& seqs, const int open, const int extend, const int match, const int mismatch) {
    // get CPU run time
    //std::cout << "running cpu" << std::endl;
    ScoreTime cpu_st = alignCPU(seqs.A, seqs.B, open, extend, match, mismatch);
    // get GPU run time
    //std::cout << "running gpu" << std::endl;
    //ScoreTime gpu_st = alignGPU(seqs.A, seqs.B, open, extend, match, mismatch);
    // get sequential run time
    //std::cout << "running sequential" << std::endl;
    //ScoreTime nw_affine = alignNW_Affine(seqs.A, seqs.B, match, mismatch, open, extend);


    // Speedups
    //float speedup_gpu_vs_cpu = cpu_st.time / gpu_st.time;
    //float speedup_seq_vs_gpu = nw_affine.time / gpu_st.time;
    //float speedup_seq_vs_cpu = nw_affine.time / cpu_st.time;

    //std::cout << "\n=== Alignment Comparison, n = " << seqs.A.size() << " ===\n";
    //std::cout << "Times in ms, CPU:" << cpu_st.time << " GPU:" << gpu_st.time << std::endl;
    
    //std::cout << "Sequential NW:  Score = " << nw_affine.score << ", Time = " << nw_affine.time << " ms\n";
    std::cout << "Gotoh (CPU):    Score = " << cpu_st.score << ", Time = " << cpu_st.time << " ms\n";
    std::cout << "Gotoh (GPU):    Score = " << gpu_st.score << ", Time = " << gpu_st.time << " ms\n";
    
    std::cout << "\n=== Speedups ===\n";
    //std::cout << "GPU vs CPU:     " << speedup_gpu_vs_cpu << "x\n";
    //std::cout << "GPU vs NW:      " << speedup_seq_vs_gpu << "x\n";
    //std::cout << "CPU vs NW:      " << speedup_seq_vs_cpu << "x\n";
    
}

int main() {
    /*
    std::vector<Strings> seqs = {
        {"GATTACA", "GATTACA"},
        {"GATTACA", "GATTTTA"},
        {"GATTACA", "A-TTACA"},
        {"GATTACA", "GATT---"},
        {"GATTACA", "GATTACC"},
        {"ACTGACTG", "ACT-ACTG"},
        {"AAAAAAA", "AAAGAAA"},
        {"TTTTTTT", "TTTCTTT"},
        {"CGTACGT", "CGTTCGT"},
        {"GGCATGC", "GGC-TGC"}
    };
    for (int n = 100; n < 1000; n += 100) {
        std::string A(n, 'A');
        for (int i = 0; i < n; i += 4) A[i] = 'C';

        std::string B(n, 'A');
        for (int i = 0; i < n; i += 5) B[i] = 'G';

        //std::string A = fastaReader("seqA.fasta");
        //std::string B = fastaReader("seqB.fasta");
        Strings seqs = {A,B};


        int openGap = 3;
        int extendGap = 1;
        int match = 3;
        int mismatch = -1;
        testCPUandGPU(seqs, openGap, extendGap, match, mismatch);
    }
     */

    //std::cout << "Running CPU tests:\n";
    //BatchTestCPU(seqs, openGap, extendGap, match, mismatch);

    //std::cout << "\nRunning GPU tests:\n";
    //BatchTestGPU(seqs, openGap, extendGap, match, mismatch);

    // running every pair of sequences on CPU and GPU to compare runtimes
    //for (int i = 0; i < seqs.size(); i++) {
    //    testCPUandGPU(seqs[i], openGap, extendGap, match, mismatch);
    //}

    std::string file = "uniref50.fasta";
    int openGap = 5;
    int extendGap = 1;
    int match = 3;
    int mismatch = -1;
    int length = 50000;

    //for (int i = 0; i < 50; i++) {
    std::string seq1 = getFastaSequenceByIndex(file, 1);
    std::string seq2 = getFastaSequenceByIndex(file, 2); 

    string seqA = seq1.substr(0, length); 
    string seqB = seq2.substr(0, length);
    //std::cout << seq1 << std::endl;
    //std::cout << seq2 << std::endl;
    Strings seqs = {seqA,seqB};
        
    testCPUandGPU(seqs, openGap, extendGap, match, mismatch);
    //}

    return 0;
}

// Example output: Scores (CPU/GPU): 208  208 Times taken: 14.2469 ms (CPU) 8.09587 ms (GPU). Speedup: 1.75978