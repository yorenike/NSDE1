#ifndef FD_DISCRETIZATION_H
#define FD_DISCRETIZATION_H

#include <vector>
#include <string>
#include "Grid.h"
#include "TestFunction.h"

// 边界条件结构体
struct BoundaryConfig {
    std::string type;      // "dirichlet" 或 "neumann"
    std::string value;     // 表达式或 "from_test_function"
    bool from_test;        // 是否从测试函数获取
};

class FDDiscretization {
private:
    const Grid& grid;
    const TestFunction& test_func;
    
    // 四个边界的配置
    BoundaryConfig bc_left;
    BoundaryConfig bc_right;
    BoundaryConfig bc_bottom;
    BoundaryConfig bc_top;
    
    // 线性系统
    std::vector<std::vector<double>> A;
    std::vector<double> F;
    std::vector<double> U;
    
    int N;
    
    void initializeSystem();
    
    // 边界条件判断
    bool isNeumann(const std::string& side) const;
    bool isDirichlet(const std::string& side) const;
    bool isDirichletOnPoint(int i, int j) const;
    bool isNeumannOnPoint(int i, int j) const;
    bool isCornerPoint(int i, int j) const;
    
    // 获取边界值
    double getBoundaryValue(int i, int j, const std::string& side) const;
    double getNeumannValue(int i, int j, const std::string& side) const;
    
    // 组装各种类型的点
    void assembleDirichletPoint(int i, int j, int idx);
    void assembleNeumannBoundaryPoint(int i, int j, int idx);
    void assembleNeumannCornerPoint(int i, int j, int idx);
    void assembleRegularPoint(int i, int j, int idx);
    void assembleIrregularPoint(int i, int j, int idx);
    
    // 辅助函数
    double findDistanceToHoleBoundary(int i, int j, const std::string& direction);
    
public:
    // 构造函数
    FDDiscretization(const Grid& grid, const TestFunction& test_func,
                     const BoundaryConfig& left, const BoundaryConfig& right,
                     const BoundaryConfig& bottom, const BoundaryConfig& top);
    
    void assemble();
    void solve();
    
    const std::vector<double>& getSolution() const { return U; }
    const std::vector<std::vector<double>>& getMatrix() const { return A; }
    const std::vector<double>& getRHS() const { return F; }
    
    void printSystem() const;
    void computeError(double& L1, double& L2, double& Linf) const;
    void setSolution(const std::vector<double>& solution) { U = solution; }
};

#endif // FD_DISCRETIZATION_H