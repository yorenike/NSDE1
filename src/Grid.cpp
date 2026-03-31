#include "Grid.h"
#include <iostream>
#include <cmath>
#include <algorithm>

// 构造函数：正方形区域
Grid::Grid(int nx_, int ny_)
    : nx(nx_), ny(ny_), N(0), domain_type("square"), 
      hole_cx(0.0), hole_cy(0.0), hole_r(0.0) {
    
    // 计算步长
    hx = 1.0 / (nx + 1);
    hy = 1.0 / (ny + 1);
    
    // 生成网格点
    generatePoints();
    
    // 分配全局编号
    assignGlobalIndices();
}

// 构造函数：带圆孔的正方形区域
Grid::Grid(int nx_, int ny_, double hole_cx_, double hole_cy_, double hole_r_)
    : nx(nx_), ny(ny_), N(0), domain_type("square_with_hole"),
      hole_cx(hole_cx_), hole_cy(hole_cy_), hole_r(hole_r_) {
    
    // 计算步长
    hx = 1.0 / (nx + 1);
    hy = 1.0 / (ny + 1);
    
    // 生成网格点
    generatePoints();
    
    // 分配全局编号
    assignGlobalIndices();
}

// 判断点是否在圆孔内
bool Grid::isInsideHole(double x, double y) const {
    if (domain_type != "square_with_hole") return false;
    
    double dx = x - hole_cx;
    double dy = y - hole_cy;
    return (dx * dx + dy * dy) < (hole_r * hole_r);
}

// 判断是否是边界点
// i: y方向索引 (0 到 ny+1), j: x方向索引 (0 到 nx+1)
bool Grid::isBoundaryPoint(int i, int j) const {
    return (i == 0 || i == ny + 1 || j == 0 || j == nx + 1);
}

// 判断是否是规则点
// 规则点：所有四个邻居都存在且不在圆孔内
bool Grid::isRegularPoint(int i, int j) const {
    // 必须是内部点
    if (isBoundaryPoint(i, j)) return false;
    
    // 检查点本身是否在圆孔内
    if (isInsideHole(points[i][j].x, points[i][j].y)) return false;
    
    // 检查四个邻居是否都在圆孔外
    bool south_ok = !isInsideHole(points[i-1][j].x, points[i-1][j].y);
    bool north_ok = !isInsideHole(points[i+1][j].x, points[i+1][j].y);
    bool west_ok = !isInsideHole(points[i][j-1].x, points[i][j-1].y);
    bool east_ok = !isInsideHole(points[i][j+1].x, points[i][j+1].y);
    
    return (south_ok && north_ok && west_ok && east_ok);
}

// 判断是否是不规则点
// 不规则点：内部点，不在圆孔内，但至少有一个邻居在圆孔内
bool Grid::isIrregularPoint(int i, int j) const {
    // 必须是内部点
    if (isBoundaryPoint(i, j)) return false;
    
    // 点本身不能在圆孔内
    if (isInsideHole(points[i][j].x, points[i][j].y)) return false;
    
    // 检查是否至少有一个邻居在圆孔内
    bool south_hole = isInsideHole(points[i-1][j].x, points[i-1][j].y);
    bool north_hole = isInsideHole(points[i+1][j].x, points[i+1][j].y);
    bool west_hole = isInsideHole(points[i][j-1].x, points[i][j-1].y);
    bool east_hole = isInsideHole(points[i][j+1].x, points[i][j+1].y);
    
    return (south_hole || north_hole || west_hole || east_hole);
}

// 生成网格点
// points[i][j] 中 i 是 y 方向索引（对应 y 坐标），j 是 x 方向索引（对应 x 坐标）
// 生成网格点（修正版）
void Grid::generatePoints() {
    // 调整大小：ny+2 行（y方向），每行 nx+2 列（x方向）
    points.resize(ny + 2, std::vector<GridPoint>(nx + 2));
    idx_map.resize(ny + 2, std::vector<int>(nx + 2, -1));
    
    // ========== 阶段1：设置坐标和基本类型 ==========
    for (int i = 0; i <= ny + 1; i++) {
        for (int j = 0; j <= nx + 1; j++) {
            double x = j * hx;
            double y = i * hy;
            
            points[i][j].x = x;
            points[i][j].y = y;
            
            // 判断点类型（初步）
            if (isBoundaryPoint(i, j)) {
                points[i][j].type = PointType::BOUNDARY;
            } else if (isInsideHole(x, y)) {
                points[i][j].type = PointType::HOLE;
            } else {
                // 其他点先标记为 REGULAR
                points[i][j].type = PointType::REGULAR;
            }
        }
    }
    
    // ========== 阶段2：识别不规则点 ==========
    // 不规则点：内部点，不在圆孔内，但至少有一个邻居在圆孔内
    for (int i = 1; i <= ny; i++) {
        for (int j = 1; j <= nx; j++) {
            // 跳过边界点和圆孔内点
            if (isBoundaryPoint(i, j)) continue;
            if (points[i][j].type == PointType::HOLE) continue;
            
            // 检查是否有邻居是圆孔内点
            bool has_hole_neighbor = false;
            
            if (points[i-1][j].type == PointType::HOLE) has_hole_neighbor = true;
            if (points[i+1][j].type == PointType::HOLE) has_hole_neighbor = true;
            if (points[i][j-1].type == PointType::HOLE) has_hole_neighbor = true;
            if (points[i][j+1].type == PointType::HOLE) has_hole_neighbor = true;
            
            if (has_hole_neighbor) {
                points[i][j].type = PointType::IRREGULAR;
            }
            // 否则保持为 REGULAR
        }
    }
}
// 分配全局编号
// 编号顺序：先固定 y，变化 x（即先遍历 x 方向，再遍历 y 方向）
// 这样编号与物理直觉一致：相邻的 x 方向点编号连续
void Grid::assignGlobalIndices() {
    N = 0;
    for (int i = 1; i <= ny; i++) {        // i: y方向索引（外层循环）
        for (int j = 1; j <= nx; j++) {    // j: x方向索引（内层循环）
            PointType type = points[i][j].type;
            if (type == PointType::REGULAR || type == PointType::IRREGULAR) {
                points[i][j].global_idx = N;
                idx_map[i][j] = N;
                N++;
            } else {
                points[i][j].global_idx = -1;
                idx_map[i][j] = -1;
            }
        }
    }
}

// 获取点信息
GridPoint Grid::getPoint(int i, int j) const {
    return points[i][j];
}

double Grid::getX(int i, int j) const {
    return points[i][j].x;
}

double Grid::getY(int i, int j) const {
    return points[i][j].y;
}

PointType Grid::getPointType(int i, int j) const {
    return points[i][j].type;
}

int Grid::getGlobalIdx(int i, int j) const {
    return idx_map[i][j];
}

// 获取所有边界点
std::vector<std::pair<int, int>> Grid::getBoundaryPoints() const {
    std::vector<std::pair<int, int>> boundary_points;
    
    for (int i = 0; i <= ny + 1; i++) {
        for (int j = 0; j <= nx + 1; j++) {
            if (isBoundaryPoint(i, j)) {
                boundary_points.push_back({i, j});
            }
        }
    }
    
    return boundary_points;
}

// 获取所有方程离散点（规则点 + 不规则点）
std::vector<std::pair<int, int>> Grid::getEquationPoints() const {
    std::vector<std::pair<int, int>> eq_points;
    
    for (int i = 1; i <= ny; i++) {
        for (int j = 1; j <= nx; j++) {
            PointType type = points[i][j].type;
            if (type == PointType::REGULAR || type == PointType::IRREGULAR) {
                eq_points.push_back({i, j});
            }
        }
    }
    
    return eq_points;
}

// 获取所有规则点
std::vector<std::pair<int, int>> Grid::getRegularPoints() const {
    std::vector<std::pair<int, int>> regular_points;
    
    for (int i = 1; i <= ny; i++) {
        for (int j = 1; j <= nx; j++) {
            if (points[i][j].type == PointType::REGULAR) {
                regular_points.push_back({i, j});
            }
        }
    }
    
    return regular_points;
}

// 获取所有不规则点
std::vector<std::pair<int, int>> Grid::getIrregularPoints() const {
    std::vector<std::pair<int, int>> irregular_points;
    
    for (int i = 1; i <= ny; i++) {
        for (int j = 1; j <= nx; j++) {
            if (points[i][j].type == PointType::IRREGULAR) {
                irregular_points.push_back({i, j});
            }
        }
    }
    
    return irregular_points;
}

// 打印网格信息
void Grid::printGridInfo() const {
    std::cout << "=== Grid Information ===" << std::endl;
    std::cout << "Domain type: " << domain_type << std::endl;
    std::cout << "Grid size: " << nx << " (x) x " << ny << " (y)" << std::endl;
    std::cout << "Step size: hx = " << hx << ", hy = " << hy << std::endl;
    
    if (domain_type == "square_with_hole") {
        std::cout << "Hole: center=(" << hole_cx << ", " << hole_cy 
                  << "), radius=" << hole_r << std::endl;
    }
    
    std::cout << "\nTotal equation points: " << N << std::endl;
    std::cout << "  Regular points: " << getRegularPoints().size() << std::endl;
    std::cout << "  Irregular points: " << getIrregularPoints().size() << std::endl;
    std::cout << "Boundary points: " << getBoundaryPoints().size() << std::endl;
    
    // 打印点类型矩阵（用于调试）
    // 打印时：外层循环 i（y方向，从上到下），内层循环 j（x方向，从左到右）
    std::cout << "\nPoint type matrix (R=Regular, I=Irregular, B=Boundary, H=Hole):" << std::endl;
    std::cout << "y ↑" << std::endl;
    for (int i = ny + 1; i >= 0; i--) {        // 从上到下打印（y从大到小）
        std::cout << "  ";
        for (int j = 0; j <= nx + 1; j++) {    // 从左到右打印（x从小到大）
            char c;
            switch (points[i][j].type) {
                case PointType::REGULAR:   c = 'R'; break;
                case PointType::IRREGULAR: c = 'I'; break;
                case PointType::BOUNDARY:  c = 'B'; break;
                case PointType::HOLE:      c = 'H'; break;
                default:                   c = '?'; break;
            }
            std::cout << c << " ";
        }
        std::cout << "  y=" << i * hy << std::endl;
    }
    std::cout << "  x→" << std::endl;
    std::cout << "  ";
    for (int j = 0; j <= nx + 1; j++) {
        std::cout << j * hx << " ";
    }
    std::cout << std::endl;
}