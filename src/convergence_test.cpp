#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "Grid.h"
#include "TestFunction.h"
#include "FDDiscretization.h"

// 使用内置高斯消元法（不用 LAPACK）
std::vector<double> solveGaussian(const std::vector<std::vector<double>>& A, 
                                   const std::vector<double>& F) {
    int n = A.size();
    std::vector<std::vector<double>> augmented(n, std::vector<double>(n + 1));
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][n] = F[i];
    }
    
    // 高斯消元
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

struct ConvergenceResult {
    int nx;
    double h;
    double L1_error;
    double L2_error;
    double Linf_error;
    double L1_order;
    double L2_order;
    double Linf_order;
};

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Convergence Test for Poisson Equation" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // 测试的网格尺寸
    std::vector<int> grid_sizes = {8, 16, 32, 64};
    std::vector<ConvergenceResult> results;
    
    // 测试函数
    std::string test_names[] = {"exp_y_sin_x", "polynomial", "trigonometric"};
    
    for (const auto& test_name : test_names) {
        std::cout << "\n\n========== Test Function: " << test_name << " ==========" << std::endl;
        
        TestFunction test_func(test_name);
        std::vector<ConvergenceResult> test_results;
        
        for (size_t idx = 0; idx < grid_sizes.size(); idx++) {
            int nx = grid_sizes[idx];
            int ny = nx;
            double h = 1.0 / (nx + 1);
            
            std::cout << "\n--- Grid: " << nx << " x " << ny << " (h = " << h << ") ---" << std::endl;
            
            // 创建网格（正方形区域）
            Grid grid(nx, ny);
            std::string bc_type = "dirichlet";
            
            // 离散化
            FDDiscretization fd(grid, test_func, bc_type);
            fd.assemble();
            
            // 求解
            std::vector<double> solution;
            try {
                solution = solveGaussian(fd.getMatrix(), fd.getRHS());
            } catch (const std::exception& e) {
                std::cerr << "Error: " << e.what() << std::endl;
                continue;
            }
            
            // 计算误差
            double L1, L2, Linf;
            fd.computeError(L1, L2, Linf);
            
            std::cout << std::scientific << std::setprecision(6);
            std::cout << "L1 error:   " << L1 << std::endl;
            std::cout << "L2 error:   " << L2 << std::endl;
            std::cout << "L∞ error:   " << Linf << std::endl;
            
            ConvergenceResult result;
            result.nx = nx;
            result.h = h;
            result.L1_error = L1;
            result.L2_error = L2;
            result.Linf_error = Linf;
            
            // 计算收敛阶
            if (idx > 0) {
                double h_prev = 1.0 / (grid_sizes[idx-1] + 1);
                result.L1_order = computeConvergenceOrder(
                    test_results[idx-1].L1_error, L1, h_prev, h);
                result.L2_order = computeConvergenceOrder(
                    test_results[idx-1].L2_error, L2, h_prev, h);
                result.Linf_order = computeConvergenceOrder(
                    test_results[idx-1].Linf_error, Linf, h_prev, h);
                
                std::cout << "L1 order:   " << result.L1_order << std::endl;
                std::cout << "L2 order:   " << result.L2_order << std::endl;
                std::cout << "L∞ order:   " << result.Linf_order << std::endl;
            } else {
                result.L1_order = 0;
                result.L2_order = 0;
                result.Linf_order = 0;
            }
            
            test_results.push_back(result);
        }
        
        // 打印汇总表格
        std::cout << "\n--- Convergence Summary for " << test_name << " ---" << std::endl;
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "nx      h        L1_error    L1_order    L2_error    L2_order    Linf_error  Linf_order" << std::endl;
        std::cout << "------------------------------------------------------------------------" << std::endl;
        
        for (const auto& r : test_results) {
            std::cout << std::setw(4) << r.nx << "   "
                      << std::setw(8) << r.h << "   "
                      << std::scientific << std::setprecision(6)
                      << std::setw(10) << r.L1_error << "   "
                      << std::fixed << std::setprecision(2)
                      << std::setw(8) << r.L1_order << "   "
                      << std::scientific << std::setprecision(6)
                      << std::setw(10) << r.L2_error << "   "
                      << std::fixed << std::setprecision(2)
                      << std::setw(8) << r.L2_order << "   "
                      << std::scientific << std::setprecision(6)
                      << std::setw(10) << r.Linf_error << "   "
                      << std::fixed << std::setprecision(2)
                      << std::setw(8) << r.Linf_order << std::endl;
        }
        
        // 验证收敛阶
        double expected_order = 2.0;
        bool all_second_order = true;
        for (size_t i = 1; i < test_results.size(); i++) {
            if (std::abs(test_results[i].L2_order - expected_order) > 0.5) {
                all_second_order = false;
            }
        }
        
        if (all_second_order && test_results.size() > 1) {
            std::cout << "\n✓ Second-order convergence (O(h²)) verified!" << std::endl;
        } else {
            std::cout << "\n⚠ Convergence order may not be second-order." << std::endl;
        }
    }
    
    return 0;
}