#pragma once
#ifndef GOTOH_LINEARSPACE_H
#define GOTOH_LINEARSPACE_H

#include <string>
#include <vector>
#include <tuple>


std::tuple<std::vector<int>, std::vector<int>, std::vector<int>>
compute_linear_space_cpu(
    const std::string& A,
    const std::string& B,
    int openGap,
    int extendGap,
    const std::vector<std::vector<int>>& submat,
    bool reverse = false
);
#endif