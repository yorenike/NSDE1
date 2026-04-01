#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include "JSONParser.h"
#include "Grid.h"
#include "TestFunction.h"
#include "FDDiscretization.h"

// 从 JSON 读取边界配置
BoundaryConfig readBoundaryConfig(const JSONParser& config, const std::string& side) {
    BoundaryConfig bc;
    std::string base_key = "boundary." + side;
    
    bc.type = config.getString(base_key + ".type");
    bc.value = config.getString(base_key + ".value");
    bc.from_test = (bc.value == "from_test_function");
    
    return bc;
}

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

void saveResults(const std::string& filename, const Grid& grid, 
                 const std::vector<double>& solution, const TestFunction& test_func) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Warning: Cannot open file " << filename << std::endl;
        return;
    }
    
    file << "# x y u_computed u_exact error" << std::endl;
    file << std::fixed << std::setprecision(8);
    
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
        
        file << x << " " << y << " " << u_computed << " " << u_exact << " " << error << std::endl;
    }
    
    file.close();
    std::cout << "Results saved to " << filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.json" << std::endl;
        return 1;
    }
    
    try {
        JSONParser config(argv[1]);
        
        std::string domain_type = config.getString("problem.domain_type");
        int nx = config.getInt("grid.nx");
        int ny = config.getInt("grid.ny");
        
        // 读取四个边界的配置
        BoundaryConfig bc_left = readBoundaryConfig(config, "left");
        BoundaryConfig bc_right = readBoundaryConfig(config, "right");
        BoundaryConfig bc_bottom = readBoundaryConfig(config, "bottom");
        BoundaryConfig bc_top = readBoundaryConfig(config, "top");
        
        std::string test_name = config.getString("test_function.name");
        bool verbose = config.getBoolOrDefault("output.verbose", false);
        bool save_solution = config.getBoolOrDefault("output.save_solution", false);
        std::string output_dir = config.getStringOrDefault("output.output_dir", "./results");
        
        if (verbose) {
            std::cout << "========================================" << std::endl;
            std::cout << "Poisson Equation Solver" << std::endl;
            std::cout << "========================================" << std::endl;
            std::cout << "Domain type: " << domain_type << std::endl;
            std::cout << "Grid: " << nx << " x " << ny << std::endl;
            std::cout << "Boundary conditions:" << std::endl;
            std::cout << "  Left:   " << bc_left.type << std::endl;
            std::cout << "  Right:  " << bc_right.type << std::endl;
            std::cout << "  Bottom: " << bc_bottom.type << std::endl;
            std::cout << "  Top:    " << bc_top.type << std::endl;
            std::cout << "Test function: " << test_name << std::endl;
        }
        
        // 创建网格
        Grid* grid = nullptr;
        if (domain_type == "square") {
            grid = new Grid(nx, ny);
        } else if (domain_type == "square_with_hole") {
            double hole_cx = config.getDouble("problem.hole_center_x");
            double hole_cy = config.getDouble("problem.hole_center_y");
            double hole_r = config.getDouble("problem.hole_radius");
            grid = new Grid(nx, ny, hole_cx, hole_cy, hole_r);
        } else {
            throw std::runtime_error("Unknown domain type: " + domain_type);
        }
        
        if (verbose) {
            grid->printGridInfo();
        }
        
        // 创建测试函数
        TestFunction test_func(test_name);
        
        // 离散化
        FDDiscretization fd(*grid, test_func, bc_left, bc_right, bc_bottom, bc_top);
        fd.assemble();
        
        if (verbose) {
            std::cout << "\nAssembled linear system with " << grid->getN() << " equations" << std::endl;
        }
        
        // 求解
        std::vector<double> solution;
        try {
            solution = solveGaussian(fd.getMatrix(), fd.getRHS());
            fd.setSolution(solution);
        } catch (const std::exception& e) {
            std::cerr << "Solve failed: " << e.what() << std::endl;
            return 1;
        }
        
        double L1, L2, Linf;
        fd.computeError(L1, L2, Linf);
        
        std::cout << "\n=== Error Norms ===" << std::endl;
        std::cout << std::scientific << std::setprecision(6);
        std::cout << "L1 norm:   " << L1 << std::endl;
        std::cout << "L2 norm:   " << L2 << std::endl;
        std::cout << "L∞ norm:   " << Linf << std::endl;
        
        // 保存结果
        if (save_solution) {
            std::string filename = output_dir + "/solution_nx" + std::to_string(nx) + ".dat";
            saveResults(filename, *grid, solution, test_func);
        }
        
        delete grid;
        
        std::cout << "\n✓ Solver completed successfully!" << std::endl;






       std::cout << "\n=== My Test ===" << std::endl;
       fd.printSystem();




        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}