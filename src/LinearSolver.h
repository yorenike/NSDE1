#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <vector>

namespace LinearSolver {
    // 高斯消元法
    std::vector<double> gaussianElimination(
        const std::vector<std::vector<double>>& A,
        const std::vector<double>& F);
}

#endif