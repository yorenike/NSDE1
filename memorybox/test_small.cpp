#include <iostream>
#include <iomanip>
#include "Grid.h"
#include "TestFunction.h"
#include "FDDiscretization.h"

// 高斯消元法
std::vector<double> solveGaussian(const std::vector<std::vector<double>>& A, 
                                   const std::vector<double>& F) {
    int n = A.size();
    if (n == 0) return std::vector<double>();
    
    std::vector<std::vector<double>> augmented(n, std::vector<double>(n + 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][n] = F[i];
    }
    
    for (int k = 0; k < n; k++) {
        int max_row = k;
        for (int i = k + 1; i < n; i++) {
            if (std::abs(augmented[i][k]) > std::abs(augmented[max_row][k])) {
                max_row = i;
            }
        }
        
        if (std::abs(augmented[max_row][k]) < 1e-12) {
            throw std::runtime_error("Singular matrix");
        }
        
        std::swap(augmented[k], augmented[max_row]);
        
        for (int i = k + 1; i < n; i++) {
            double factor = augmented[i][k] / augmented[k][k];
            for (int j = k; j <= n; j++) {
                augmented[i][j] -= factor * augmented[k][j];
            }
        }
    }
    
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

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Debug Test with 3x3 Grid" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // 小网格
    int nx = 3, ny = 3;
    Grid grid(nx, ny);
    double h = grid.getHx();
    
    std::cout << "\nh = " << h << ", 1/h^2 = " << 1.0/(h*h) << std::endl;
    
    // 测试函数
    TestFunction test_func("exp_y_sin_x");
    
    // 边界条件
    BoundaryConfig bc_dirichlet;
    bc_dirichlet.type = "dirichlet";
    bc_dirichlet.value = "from_test_function";
    bc_dirichlet.from_test = true;
    
    // 创建离散化
    FDDiscretization fd(grid, test_func, bc_dirichlet, bc_dirichlet, bc_dirichlet, bc_dirichlet);
    fd.assemble();
    
    // 打印矩阵和右端项
    const auto& A = fd.getMatrix();
    const auto& F = fd.getRHS();
    int N = A.size();
    
    std::cout << "\n--- Matrix A (first 5 rows) ---" << std::endl;
    for (int i = 0; i < std::min(5, N); i++) {
        for (int j = 0; j < std::min(5, N); j++) {
            std::cout << std::setw(10) << std::fixed << std::setprecision(2) << A[i][j] << " ";
        }
        std::cout << " | " << std::setw(10) << F[i] << std::endl;
    }
    
    // 求解
    std::vector<double> solution = solveGaussian(A, F);
    
    // ========== 关键：把解设置回 fd 对象 ==========
    fd.setSolution(solution);
    
    // 打印结果
    std::cout << "\n--- Solution and Exact Values ---" << std::endl;
    auto eq_points = grid.getEquationPoints();
    for (const auto& point : eq_points) {
        int i = point.first;
        int j = point.second;
        int idx = grid.getGlobalIdx(i, j);
        
        double x = grid.getX(i, j);
        double y = grid.getY(i, j);
        double u_computed = solution[idx];
        double u_exact = test_func.exact(x, y);
        double error = std::abs(u_computed - u_exact);
        
        std::cout << "(" << i << "," << j << ") x=" << x << " y=" << y 
                  << " computed=" << std::fixed << std::setprecision(6) << u_computed
                  << " exact=" << u_exact
                  << " error=" << error << std::endl;
    }
    
    // 现在 computeError 会使用正确的解
    double L1, L2, Linf;
    fd.computeError(L1, L2, Linf);
    
    std::cout << "\n--- Error Norms ---" << std::endl;
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "L1 norm:   " << L1 << std::endl;
    std::cout << "L2 norm:   " << L2 << std::endl;
    std::cout << "L∞ norm:   " << Linf << std::endl;
    
    // 验证
    if (Linf < 0.01) {
        std::cout << "\n✓ Results look reasonable! (max error = " << Linf << ")" << std::endl;
    } else {
        std::cout << "\n✗ Errors are too large! Something is wrong." << std::endl;
    }
    
    return 0;
}