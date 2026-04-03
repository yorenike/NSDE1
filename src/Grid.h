#ifndef GRID_H
#define GRID_H

#include <vector>
#include <string>
#include <functional>

enum class PointType {
    REGULAR,      // 规则点：五点差分可用
    IRREGULAR,    // 不规则点：靠近圆孔，需要特殊处理
    BOUNDARY,     // 边界点：Dirichlet条件直接给定
    HOLE          // 圆孔内的点：不参与计算
};

struct GridPoint {
    double x;           // x 坐标
    double y;           // y 坐标
    PointType type;     // 点类型
    int global_idx;     // 全局编号（-1 表示不是方程离散点）
    
    GridPoint() : x(0.0), y(0.0), type(PointType::BOUNDARY), global_idx(-1) {}
    GridPoint(double x_, double y_) : x(x_), y(y_), type(PointType::BOUNDARY), global_idx(-1) {}
};

class Grid {
private:
    int nx;                     // x 方向内部点数（不含边界）
    int ny;                     // y 方向内部点数（不含边界）
    int N;                      // 总方程离散点数（所有非HOLE点）
    double hx;                  // x 方向步长
    double hy;                  // y 方向步长
    
    // 区域参数
    std::string domain_type;    // "square" 或 "square_with_hole"
    double hole_cx, hole_cy;    // 圆孔中心坐标
    double hole_r;              // 圆孔半径
    
    // 网格数据：points[i][j] 中 i 是 y 方向索引，j 是 x 方向索引
    std::vector<std::vector<GridPoint>> points;  // points[i][j] 对应 (j*hx, i*hy)
    std::vector<std::vector<int>> idx_map;       // -1 表示不是方程离散点
    
    // 辅助函数
    bool validateHoleParameters() const;
    bool isInsideHole(double x, double y) const;
    bool isBoundaryPoint(int i, int j) const;    // i: y索引, j: x索引
    bool isRegularPoint(int i, int j) const;
    bool isIrregularPoint(int i, int j) const;
    void generatePoints();
    void assignGlobalIndices();
    
public:
    // 构造函数：正方形区域
    Grid(int nx_, int ny_);
    
    // 构造函数：带圆孔的正方形区域
    Grid(int nx_, int ny_, double hole_cx_, double hole_cy_, double hole_r_);
    
    // 获取网格参数
    int getNx() const { return nx; }
    int getNy() const { return ny; }
    int getN() const { return N; }
    double getHx() const { return hx; }
    double getHy() const { return hy; }
    
    // 获取圆孔参数
    bool hasHole() const { return domain_type == "square_with_hole"; }
    double getHoleCx() const { return hole_cx; }
    double getHoleCy() const { return hole_cy; }
    double getHoleR() const { return hole_r; }
    
    // 获取点信息：i 是 y 方向索引（0 到 ny+1），j 是 x 方向索引（0 到 nx+1）
    GridPoint getPoint(int i, int j) const;
    double getX(int i, int j) const;  // i: y索引, j: x索引
    double getY(int i, int j) const;
    PointType getPointType(int i, int j) const;
    int getGlobalIdx(int i, int j) const;
    
    // 获取边界点（用于 Dirichlet 条件）
    std::vector<std::pair<int, int>> getBoundaryPoints() const;
    
    // 获取所有方程离散点（所有非HOLE点）
    std::vector<std::pair<int, int>> getEquationPoints() const;
    
    // 获取规则点
    std::vector<std::pair<int, int>> getRegularPoints() const;
    
    // 获取不规则点
    std::vector<std::pair<int, int>> getIrregularPoints() const;
    
    // 调试：打印网格信息
    void printGridInfo() const;
};

#endif // GRID_H