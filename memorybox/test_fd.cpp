#include <iostream>
#include <iomanip>
#include "Grid.h"
#include "TestFunction.h"
#include "FDDiscretization.h"

int main() {
    std::cout << "========================================" << std::endl;
    std::cout << "Test FDDiscretization with 3x3 Grid" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // ========== 1. 创建网格（3x3 正方形区域）==========
    int nx = 3;
    int ny = 3;
    Grid grid(nx, ny);
    
    double h = grid.getHx();
    double inv_h2 = 1.0 / (h * h);
    
    std::cout << "\nGrid parameters:" << std::endl;
    std::cout << "  nx = " << nx << ", ny = " << ny << std::endl;
    std::cout << "  h = " << h << std::endl;
    std::cout << "  1/h^2 = " << inv_h2 << std::endl;
    std::cout << "  4/h^2 = " << 4.0 * inv_h2 << std::endl;
    std::cout << "  -1/h^2 = " << -inv_h2 << std::endl;
    
    // ========== 2. 创建测试函数 ==========
    TestFunction test_func("polynomial");
    
    // ========== 3. 创建 FDDiscretization 并组装 ==========
    std::string bc_type = "dirichlet";
    FDDiscretization fd(grid, test_func, bc_type);
    
    std::cout << "\n--- Assembling linear system ---" << std::endl;
    fd.assemble();
    
    // ========== 4. 获取矩阵并打印 ==========
    const auto& A = fd.getMatrix();
    const auto& F = fd.getRHS();
    int N = A.size();
    
    std::cout << "\n--- Matrix A (" << N << "x" << N << ") ---" << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << std::setw(8) << A[i][j] << " ";
        }
        std::cout << " | " << std::setw(8) << F[i] << std::endl;
    }
    
    // ========== 5. 输出网格点索引映射 ==========
    std::cout << "\n--- Global Index Mapping (i=y-index, j=x-index) ---" << std::endl;
    for (int i = 1; i <= ny; i++) {
        for (int j = 1; j <= nx; j++) {
            int idx = grid.getGlobalIdx(i, j);
            double x = grid.getX(i, j);
            double y = grid.getY(i, j);
            std::cout << "  (" << i << "," << j << ") -> idx=" << idx 
                      << ", (x=" << std::fixed << std::setprecision(3) << x 
                      << ", y=" << y << ")" << std::endl;
        }
    }
    
    // ========== 6. 验证矩阵结构 ==========
    std::cout << "\n--- Verification ---" << std::endl;
    
    double expected_diag = 4.0 * inv_h2;
    double expected_off = -inv_h2;
    
    std::cout << "Expected diagonal: " << expected_diag << std::endl;
    std::cout << "Expected off-diagonal (neighbors): " << expected_off << std::endl;
    
    // 定义 3x3 网格的邻居关系
    // 编号：(1,1)=0, (1,2)=1, (1,3)=2
    //       (2,1)=3, (2,2)=4, (2,3)=5
    //       (3,1)=6, (3,2)=7, (3,3)=8
    int neighbors[9][4] = {
        {1, 3, -1, -1},     // 0: 右(1), 下(3)
        {0, 2, 4, -1},      // 1: 左(0), 右(2), 下(4)
        {1, 5, -1, -1},     // 2: 左(1), 下(5)
        {0, 4, 6, -1},      // 3: 上(0), 右(4), 下(6)
        {1, 3, 5, 7},       // 4: 上(1), 左(3), 右(5), 下(7)
        {2, 4, 8, -1},      // 5: 上(2), 左(4), 下(8)
        {3, 7, -1, -1},     // 6: 上(3), 右(7)
        {4, 6, 8, -1},      // 7: 上(4), 左(6), 右(8)
        {5, 7, -1, -1}      // 8: 上(5), 左(7)
    };
    
    bool correct = true;
    
    // 检查对角线
    for (int i = 0; i < N; i++) {
        if (std::abs(A[i][i] - expected_diag) > 1e-8) {
            std::cout << "ERROR: A[" << i << "][" << i << "] = " << A[i][i] 
                      << ", expected " << expected_diag << std::endl;
            correct = false;
        }
    }
    
    // 检查邻居
    for (int i = 0; i < N; i++) {
        for (int k = 0; k < 4; k++) {
            int j = neighbors[i][k];
            if (j >= 0) {
                if (std::abs(A[i][j] - expected_off) > 1e-8) {
                    std::cout << "ERROR: A[" << i << "][" << j << "] = " << A[i][j]
                              << ", expected " << expected_off << std::endl;
                    correct = false;
                }
                if (std::abs(A[i][j] - A[j][i]) > 1e-8) {
                    std::cout << "ERROR: Matrix not symmetric at (" << i << "," << j << ")" << std::endl;
                    correct = false;
                }
            }
        }
    }
    
    // 检查非邻居应该是 0
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i != j) {
                bool is_neighbor = false;
                for (int k = 0; k < 4; k++) {
                    if (neighbors[i][k] == j) {
                        is_neighbor = true;
                        break;
                    }
                }
                if (!is_neighbor && std::abs(A[i][j]) > 1e-8) {
                    std::cout << "ERROR: A[" << i << "][" << j << "] = " << A[i][j]
                              << " should be 0 (not neighbors)" << std::endl;
                    correct = false;
                }
            }
        }
    }
    
    if (correct) {
        std::cout << "\n✓ Matrix structure is CORRECT!" << std::endl;
        std::cout << "  - Diagonal: 4/h^2 = " << expected_diag << std::endl;
        std::cout << "  - Neighbors: -1/h^2 = " << expected_off << std::endl;
        std::cout << "  - Matrix is symmetric" << std::endl;
    } else {
        std::cout << "\n✗ Matrix structure has ERRORS!" << std::endl;
    }
    
    return 0;
}