#ifndef GOTOHCPU_H
#define GOTOHCPU_H
#include "types.h"
#include <string>

// necessary header file to be able to call CUDA and CPU implementations in one function call

ScoreTime alignCPU(const std::string& A, const std::string& B,
    const int openGap, const int extendGap,
    const int match, const int mismatch, int n_threads);

#endif // GOTOHCPU_H