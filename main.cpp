#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "gotohCPU.h"
#include "gotohCUDA.h"
#include "types.h"
#include<chrono>
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

/*


The score return types aren't correct here
I'll edit and revisit these functions when I want to test like 100s of sequences at once, for now I'm comparing one sequence pair at a time

void BatchTestCPU(const std::vector<Strings>& seqs, const int open, const int extend, const int match, const int mismatch) {
    for (size_t i = 0; i < seqs.size(); ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        int score = alignCPU(seqs[i].A, seqs[i].B, open, extend, match, mismatch);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;

        std::cout << "test:  " << i+1 << " score (CPU): " << score 
                  << " time taken: " << elapsed.count() << " ms" << std::endl;
    }
}

void BatchTestGPU(const std::vector<Strings>& seqs, const int open, const int extend, const int match, const int mismatch) {
    for (size_t i = 0; i < seqs.size(); ++i) {
        auto start = std::chrono::high_resolution_clock::now();
        int score = alignGPU(seqs[i].A, seqs[i].B, open, extend, match, mismatch);

        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double, std::milli> elapsed = end - start;

        std::cout << "test: " << i+1 << " score (GPU): " << score 
                  << " time taken: " << elapsed.count() << " ms" << std::endl;
    }
}
*/

// compare the CPU and GPU run times for one pair of sequences
void testCPUandGPU(const Strings& seqs, const int open, const int extend, const int match, const int mismatch) {
        // get CPU run time
        ScoreTime cpu_st = alignCPU(seqs.A, seqs.B, open, extend, match, mismatch);
        // get GPU run time
        ScoreTime gpu_st = alignGPU(seqs.A, seqs.B, open, extend, match, mismatch);

        float speedup = cpu_st.time/ gpu_st.time;

        std::cout <<"Scores (CPU/GPU): " << cpu_st.score << "  " << gpu_st.score
                << " Times taken: " << cpu_st.time << " ms (CPU) " << gpu_st.time << " ms (GPU). Speedup: " << speedup << std::endl;

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
    };*/

    std::string A = fastaReader("seqA.fasta");
    std::string B = fastaReader("seqB.fasta");
    std::vector<Strings> seqs = {{A,B}};

    int openGap = 10;
    int extendGap = 1;
    int match = 3;
    int mismatch = -1;

    //std::cout << "Running CPU tests:\n";
    //BatchTestCPU(seqs, openGap, extendGap, match, mismatch);

    //std::cout << "\nRunning GPU tests:\n";
    //BatchTestGPU(seqs, openGap, extendGap, match, mismatch);

    // running every pair of sequences on CPU and GPU to compare runtimes
    for (int i = 0; i < seqs.size(); i++) {
        testCPUandGPU(seqs[i], openGap, extendGap, match, mismatch);
    }

    return 0;
}