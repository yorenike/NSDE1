#include "LinearSolver.h"
#include <cmath>
#include <vector>
#include <stdexcept>
#include <algorithm>

namespace LinearSolver {

std::vector<double> gaussianElimination(
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& F) {
    
    int n = A.size();
    if (n == 0) return std::vector<double>();
    
    // 检查 A 和 F 的维度是否匹配
    if (F.size() != (size_t)n) {
        throw std::runtime_error("Matrix and RHS dimension mismatch");
    }
    
    // 构建增广矩阵
    std::vector<std::vector<double>> augmented(n, std::vector<double>(n + 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][n] = F[i];
    }
    
    // 高斯消元（选主元）
    for (int k = 0; k < n; k++) {
        // 选主元：找到当前列绝对值最大的行
        int max_row = k;
        for (int i = k + 1; i < n; i++) {
            if (std::abs(augmented[i][k]) > std::abs(augmented[max_row][k])) {
                max_row = i;
            }
        }
        
        // 检查主元是否为零（矩阵奇异）
        if (std::abs(augmented[max_row][k]) < 1e-12) {
            throw std::runtime_error("Singular matrix");
        }
        
        // 交换当前行与主元行
        if (max_row != k) {
            std::swap(augmented[k], augmented[max_row]);
        }
        
        // 消元
        for (int i = k + 1; i < n; i++) {
            double factor = augmented[i][k] / augmented[k][k];
            for (int j = k; j <= n; j++) {
                augmented[i][j] -= factor * augmented[k][j];
            }
        }
    }
    
    // 回代求解
    std::vector<double> U(n);
    for (int i = n - 1; i >= 0; i--) {
        U[i] = augmented[i][n];
        for (int j = i + 1; j < n; j++) {
            U[i] -= augmented[i][j] * U[j];
        }
        U[i] /= augmented[i][i];
    }
    
    return U;
}

} // namespace LinearSolver