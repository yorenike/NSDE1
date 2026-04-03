#include <iostream>
#include "Grid.h"

int main() {
    // 测试正方形区域
    std::cout << "=== Test 1: Square Domain ===" << std::endl;
    Grid square_grid(8, 8);
    square_grid.printGridInfo();
    
    std::cout << "\n";
    
    // 测试带圆孔的正方形区域
    std::cout << "=== Test 2: Square with Hole ===" << std::endl;
    Grid hole_grid(16, 16, 0.5, 0.5, 0.15);
    hole_grid.printGridInfo();
    
    // 测试获取方程点
    std::cout << "\n=== Test 3: Equation Points ===" << std::endl;
    auto eq_points = hole_grid.getEquationPoints();
    std::cout << "First 10 equation points:" << std::endl;
    for (size_t k = 0; k < std::min(size_t(20), eq_points.size()); k++) {
        int i = eq_points[k].first;
        int j = eq_points[k].second;
        int idx = hole_grid.getGlobalIdx(i, j);
        std::cout << "  (" << i << "," << j << ") -> global idx=" << idx 
                  << ", x=" << hole_grid.getX(i, j) 
                  << ", y=" << hole_grid.getY(i, j) << std::endl;
    }
    
    return 0;
}