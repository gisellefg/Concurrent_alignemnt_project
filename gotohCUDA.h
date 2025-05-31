#ifndef GOTOHCUDA_H
#define GOTOHCUDA_H
#include "types.h"
#include <string>

// necessary header file to be able to call CUDA and CPU implementations in one function call

struct ScoreTime {
    int score = 0;
    float time = 0;
    ScoreTime() = default;
    ScoreTime(int s, double t) : score(s), time(t) {}
}

ScoreTime alignGPU(const std::string& A, const std::string& B, const int openGap, const int extendGap, const int match, const int mismatch);


#endif // GOTOHCUDA_H