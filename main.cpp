#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <thread>
#include "gotohCPU.h"
#include "gotohCUDA.h"
#include "types.h"
#include <chrono>
#include "fasta_reader.cpp"
#include <iostream>
#include <random>
#include <unordered_set>
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
/*
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
*/
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



// compare the CPU and GPU run times for ONE pair of sequences (not a batch test)
void testCPUandGPU(const Strings& seqs, const int open, const int extend, const int match, const int mismatch, int n_threads) {
        // get CPU run time
        //std::cout << "running cpu" << std::endl;
        ScoreTime cpu_st = alignCPU(seqs.A, seqs.B, open, extend, match, mismatch, n_threads);
        // get GPU run time
        //std::cout << "running gpu" << std::endl;
        ScoreTime gpu_st = alignGPU(seqs.A, seqs.B, open, extend, match, mismatch);

        // Speedups
        float speedup_gpu_vs_cpu = cpu_st.time / gpu_st.time;
        std::cout << "\n=== Alignment Comparison, n = " << seqs.A.size() << " ===\n";

        std::cout << "Gotoh (CPU):    Score = " << cpu_st.score << ", Time = " << cpu_st.time << " ms\n";
        std::cout << "Gotoh (GPU):    Score = " << gpu_st.score << ", Time = " << gpu_st.time << " ms\n";
        
        std::cout << "\n=== Speedups ===\n";
        std::cout << "GPU vs CPU:     " << speedup_gpu_vs_cpu << "x\n";
}

// compare the CPU and GPU run times for ONE pair of sequences (not a batch test)
double testThreadCount(const Strings& seqs, const int open, const int extend, const int match, const int mismatch, int n_threads) {
        // get CPU run time
        //std::cout << "running cpu" << std::endl;
        ScoreTime res = alignCPU(seqs.A, seqs.B, open, extend, match, mismatch, n_threads);
        double plop = res.time;
        return plop;
        //std::cout << cpu_st.time << std::endl;
}

double testKernel(const Strings& seqs, const int open, const int extend, const int match, const int mismatch) {
    ScoreTime res = alignGPU(seqs.A, seqs.B, open, extend, match, mismatch);
    double plop = res.time;
    return plop;

}

// Code generated by ChatGPT to study how input similarity affects algorithm performance
// the code takes as input a sequence, how many characters you want to (randomly) change, and a string of characters
// to choose from
std::string mutate_amino_acids(const std::string& sequence, int num_mutations, const std::string& amino_acids) {
    if (sequence.empty() || num_mutations <= 0) {
        return sequence;
    }

    std::string mutated = sequence;
    int len = mutated.size();
    num_mutations = std::min(num_mutations, len);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::unordered_set<int> mutated_indices;

    while (mutated_indices.size() < static_cast<size_t>(num_mutations)) {
        mutated_indices.insert(gen() % len);
    }

    for (int idx : mutated_indices) {
        char original = mutated[idx];

        // Create a pool excluding the original character
        std::vector<char> options;
        for (char c : amino_acids) {
            if (c != original) {
                options.push_back(c);
            }
        }

        if (!options.empty()) {
            std::uniform_int_distribution<> dis(0, options.size() - 1);
            mutated[idx] = options[dis(gen)];
        }
    }

    return mutated;
}
int main() {

    // getting database
    std::string file = "uniref50.fasta";

    //defining penalties
    int openGap = 5;
    int extendGap = 1;
    int match = 3;
    int mismatch = -1;

    // length of sequence to study
    int length = 5000;

    // single sequence runtime comparison for all algorithms:

    // get protein sequences from database
    std::string seq1 = getFastaSequenceByIndex(file, 1); // can choose index here to be any, defined which sequence from the file it reads
    std::string seq2 = getFastaSequenceByIndex(file, 2);

    // cut to length
    string seqA = seq1.substr(0, length); 
    string seqB = seq2.substr(0, length); 
    Strings seqs = {seqA,seqB};

    // run test
    testCPUandGPU(seqs, openGap, extendGap, match, mismatch, 20); // last number here is thread count



    // IMPLEMENTATION FOR RUNNING 50 RUNS AT ONCE / REPLACING CHARACTERS
    // below

    /*

    std::string amino_acids = "ARNDCQEGHILKMFPSTWYVBZXOUJ";

    std::cout << "gpu count" << std::endl;
    for (int i = 0; i < 50; i++) {
        if (i % 10 == 0) {
        std::cout << "iteration " << i << std::endl; };

            std::string seq1 = getFastaSequenceByIndex(file, i);
            //std::string seq2 = getFastaSequenceByIndex(file, 2*i);
            //std::string seq2 = mutate_amino_acids(seq1, 0, amino_acids);

            //string seqA = seq1.substr(0, length); 
            string seqA = seq1.substr(0, length); 
            string seqB = mutate_amino_acids(seqA, change, amino_acids);
            //std::cout << seqA << " vs " << seqB << std::endl;
            
            //string seqB = seq2.substr(0, length);
            //std::cout << seq1 << std::endl;
            //std::cout << seq2 << std::endl;
            Strings seqs = {seqA,seqB};
                
            double add = testKernel(seqs, openGap, extendGap, match, mismatch);
            //testCPUandGPU(seqs, openGap, extendGap, match, mismatch, 1);
            sum += add;
        }
    sum = sum/50;
    std::cout << sum << std::endl;
    

    double sum1 = 0;
    std::cout << "cpu count" << std::endl;
    for (int i = 0; i < 50; i++) {
        if (i % 10 == 0) {
        std::cout << "iteration " << i << std::endl; };
            std::string seq1 = getFastaSequenceByIndex(file, i);
            //std::string seq2 = getFastaSequenceByIndex(file, 2*i);
            string seqA = seq1.substr(0, length); 
            string seqB = mutate_amino_acids(seqA, change, amino_acids);
           // string seqA = seq1.substr(0, length); 
            //string seqB = seq2.substr(0, length);
            //std::cout << seq1 << std::endl;
            //std::cout << seq2 << std::endl;
            Strings seqs = {seqA,seqB};
                
            double add = testThreadCount(seqs, openGap, extendGap, match, mismatch,20);
            //testCPUandGPU(seqs, openGap, extendGap, match, mismatch, 1);
            sum1 += add;
        }
    sum1 = sum1/50;
    std::cout << sum1 << std::endl;


    double sum2 = 0;
    std::cout << "sequential count" << std::endl;
    for (int i = 0; i < 50; i++) {
        if (i % 10 == 0) {
        std::cout << "iteration " << i << std::endl; };
            std::string seq1 = getFastaSequenceByIndex(file, i);
            //std::string seq2 = getFastaSequenceByIndex(file, 2*i);
            //std::string seq2 = mutate_amino_acids(seq1, 0, amino_acids);

            //string seqA = seq1.substr(0, length); 
            //string seqB = seq2.substr(0, length);
            string seqA = seq1.substr(0, length); 
            string seqB = mutate_amino_acids(seqA, change, amino_acids);
            //std::cout << seq1 << std::endl;
            //std::cout << seq2 << std::endl;
            Strings seqs = {seqA,seqB};
                
            double add = testThreadCount(seqs, openGap, extendGap, match, mismatch,1);
            //testCPUandGPU(seqs, openGap, extendGap, match, mismatch, 1);
            sum2 += add;
        }
    sum2 = sum2/50;
    std::cout << sum2 << std::endl;
    */

    return 0;
}

// Example output: Scores (CPU/GPU): 208  208 Times taken: 14.2469 ms (CPU) 8.09587 ms (GPU). Speedup: 1.75978