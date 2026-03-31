#include "FDDiscretization.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>

FDDiscretization::FDDiscretization(const Grid& grid_, const TestFunction& test_func_,
                                   const BoundaryConfig& left_, const BoundaryConfig& right_,
                                   const BoundaryConfig& bottom_, const BoundaryConfig& top_)
    : grid(grid_), test_func(test_func_),
      bc_left(left_), bc_right(right_), bc_bottom(bottom_), bc_top(top_) {
    N = grid.getN();
    initializeSystem();
}

void FDDiscretization::initializeSystem() {
    A.resize(N, std::vector<double>(N, 0.0));
    F.resize(N, 0.0);
    U.resize(N, 0.0);
}

double FDDiscretization::getBoundaryValue(int i, int j, const std::string& side) const {
    double x = grid.getX(i, j);
    double y = grid.getY(i, j);
    
    const BoundaryConfig* bc = nullptr;
    if (side == "left") bc = &bc_left;
    else if (side == "right") bc = &bc_right;
    else if (side == "bottom") bc = &bc_bottom;
    else if (side == "top") bc = &bc_top;
    else return 0.0;
    
    if (bc->from_test) {
        return test_func.getDirichletBC(x, y, side);
    } else {
        // 这里可以扩展支持表达式解析
        // 目前简单返回 0
        return 0.0;
    }
}

void FDDiscretization::assembleRegularPoint(int i, int j, int idx) {
    double h = grid.getHx();
    double h2 = h * h;
    
    double x = grid.getX(i, j);
    double y = grid.getY(i, j);
    
    A[idx][idx] = 4.0 / h2;
    
    // 西边
    if (grid.getPointType(i, j-1) == PointType::BOUNDARY) {
        // 边界点，贡献到右端项
        double bc_val = getBoundaryValue(i, j-1, "left");
        F[idx] += bc_val / h2;
    } else if (grid.getPointType(i, j-1) != PointType::HOLE) {
        int neighbor_idx = grid.getGlobalIdx(i, j-1);
        if (neighbor_idx >= 0) {
            A[idx][neighbor_idx] = -1.0 / h2;
        }
    }
    
    // 东边
    if (grid.getPointType(i, j+1) == PointType::BOUNDARY) {
        double bc_val = getBoundaryValue(i, j+1, "right");
        F[idx] += bc_val / h2;
    } else if (grid.getPointType(i, j+1) != PointType::HOLE) {
        int neighbor_idx = grid.getGlobalIdx(i, j+1);
        if (neighbor_idx >= 0) {
            A[idx][neighbor_idx] = -1.0 / h2;
        }
    }
    
    // 南边
    if (grid.getPointType(i-1, j) == PointType::BOUNDARY) {
        double bc_val = getBoundaryValue(i-1, j, "bottom");
        F[idx] += bc_val / h2;
    } else if (grid.getPointType(i-1, j) != PointType::HOLE) {
        int neighbor_idx = grid.getGlobalIdx(i-1, j);
        if (neighbor_idx >= 0) {
            A[idx][neighbor_idx] = -1.0 / h2;
        }
    }
    
    // 北边
    if (grid.getPointType(i+1, j) == PointType::BOUNDARY) {
        double bc_val = getBoundaryValue(i+1, j, "top");
        F[idx] += bc_val / h2;
    } else if (grid.getPointType(i+1, j) != PointType::HOLE) {
        int neighbor_idx = grid.getGlobalIdx(i+1, j);
        if (neighbor_idx >= 0) {
            A[idx][neighbor_idx] = -1.0 / h2;
        }
    }
    
    // 右端项：f(x,y)
    F[idx] += test_func.f(x, y);
}

void FDDiscretization::assembleIrregularPoint(int i, int j, int idx) {
    // 不规则点处理（类似规则点，但需要考虑圆孔边界）
    double h = grid.getHx();
    double h2 = h * h;
    double x = grid.getX(i, j);
    double y = grid.getY(i, j);
    
    double center_coeff = 0.0;
    double rhs = test_func.f(x, y);
    
    // 西边
    if (grid.getPointType(i, j-1) == PointType::HOLE) {
        double dist = findDistanceToHoleBoundary(i, j, "west");
        double theta = dist / h;
        double coeff = 2.0 / (theta * h2);
        center_coeff += coeff * (1.0 + theta);
        // 圆孔边界值（Dirichlet，从测试函数获取）
        double x_boundary = x - dist;
        double y_boundary = y;
        rhs += coeff * theta * test_func.exact(x_boundary, y_boundary);
    } else if (grid.getPointType(i, j-1) == PointType::BOUNDARY) {
        double bc_val = getBoundaryValue(i, j-1, "left");
        rhs += bc_val / h2;
        center_coeff += 1.0 / h2;
    } else if (grid.getPointType(i, j-1) != PointType::HOLE) {
        int neighbor_idx = grid.getGlobalIdx(i, j-1);
        if (neighbor_idx >= 0) {
            A[idx][neighbor_idx] = -1.0 / h2;
            center_coeff += 1.0 / h2;
        }
    }
    
    // 东边
    if (grid.getPointType(i, j+1) == PointType::HOLE) {
        double dist = findDistanceToHoleBoundary(i, j, "east");
        double theta = dist / h;
        double coeff = 2.0 / (theta * h2);
        center_coeff += coeff * (1.0 + theta);
        double x_boundary = x + dist;
        double y_boundary = y;
        rhs += coeff * theta * test_func.exact(x_boundary, y_boundary);
    } else if (grid.getPointType(i, j+1) == PointType::BOUNDARY) {
        double bc_val = getBoundaryValue(i, j+1, "right");
        rhs += bc_val / h2;
        center_coeff += 1.0 / h2;
    } else if (grid.getPointType(i, j+1) != PointType::HOLE) {
        int neighbor_idx = grid.getGlobalIdx(i, j+1);
        if (neighbor_idx >= 0) {
            A[idx][neighbor_idx] = -1.0 / h2;
            center_coeff += 1.0 / h2;
        }
    }
    
    // 南边
    if (grid.getPointType(i-1, j) == PointType::HOLE) {
        double dist = findDistanceToHoleBoundary(i, j, "south");
        double alpha = dist / h;
        double coeff = 2.0 / (alpha * h2);
        center_coeff += coeff * (1.0 + alpha);
        double x_boundary = x;
        double y_boundary = y - dist;
        rhs += coeff * alpha * test_func.exact(x_boundary, y_boundary);
    } else if (grid.getPointType(i-1, j) == PointType::BOUNDARY) {
        double bc_val = getBoundaryValue(i-1, j, "bottom");
        rhs += bc_val / h2;
        center_coeff += 1.0 / h2;
    } else if (grid.getPointType(i-1, j) != PointType::HOLE) {
        int neighbor_idx = grid.getGlobalIdx(i-1, j);
        if (neighbor_idx >= 0) {
            A[idx][neighbor_idx] = -1.0 / h2;
            center_coeff += 1.0 / h2;
        }
    }
    
    // 北边
    if (grid.getPointType(i+1, j) == PointType::HOLE) {
        double dist = findDistanceToHoleBoundary(i, j, "north");
        double alpha = dist / h;
        double coeff = 2.0 / (alpha * h2);
        center_coeff += coeff * (1.0 + alpha);
        double x_boundary = x;
        double y_boundary = y + dist;
        rhs += coeff * alpha * test_func.exact(x_boundary, y_boundary);
    } else if (grid.getPointType(i+1, j) == PointType::BOUNDARY) {
        double bc_val = getBoundaryValue(i+1, j, "top");
        rhs += bc_val / h2;
        center_coeff += 1.0 / h2;
    } else if (grid.getPointType(i+1, j) != PointType::HOLE) {
        int neighbor_idx = grid.getGlobalIdx(i+1, j);
        if (neighbor_idx >= 0) {
            A[idx][neighbor_idx] = -1.0 / h2;
            center_coeff += 1.0 / h2;
        }
    }
    
    A[idx][idx] = center_coeff;
    F[idx] = rhs;
}

void FDDiscretization::applyBoundaryConditions() {
    // 边界条件已在 assembleRegularPoint 和 assembleIrregularPoint 中处理
    // 这里不需要额外操作
}

void FDDiscretization::assemble() {
    initializeSystem();
    
    auto eq_points = grid.getEquationPoints();
    
    for (const auto& point : eq_points) {
        int i = point.first;
        int j = point.second;
        int idx = grid.getGlobalIdx(i, j);
        
        if (grid.getPointType(i, j) == PointType::REGULAR) {
            assembleRegularPoint(i, j, idx);
        } else if (grid.getPointType(i, j) == PointType::IRREGULAR) {
            assembleIrregularPoint(i, j, idx);
        }
    }
}

double FDDiscretization::findDistanceToHoleBoundary(int i, int j, const std::string& direction) {
    double x = grid.getX(i, j);
    double y = grid.getY(i, j);
    double h = grid.getHx();
    
    // 圆孔参数（实际应该从 grid 获取）
    double cx = 0.5, cy = 0.5, r = 0.25;
    
    double step = h / 100.0;
    double distance = 0.0;
    
    while (distance < h) {
        double test_x = x, test_y = y;
        if (direction == "west") test_x = x - distance;
        else if (direction == "east") test_x = x + distance;
        else if (direction == "south") test_y = y - distance;
        else if (direction == "north") test_y = y + distance;
        
        double dx = test_x - cx;
        double dy = test_y - cy;
        if (dx*dx + dy*dy < r*r) {
            return distance;
        }
        distance += step;
    }
    
    return h;
}

void FDDiscretization::solve() {
    if (N == 0) return;
    
    // 高斯消元法
    std::vector<std::vector<double>> augmented(N, std::vector<double>(N + 1));
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][N] = F[i];
    }
    
    for (int k = 0; k < N; k++) {
        int max_row = k;
        for (int i = k + 1; i < N; i++) {
            if (std::abs(augmented[i][k]) > std::abs(augmented[max_row][k])) {
                max_row = i;
            }
        }
        
        if (std::abs(augmented[max_row][k]) < 1e-12) {
            throw std::runtime_error("Singular matrix");
        }
        
        std::swap(augmented[k], augmented[max_row]);
        
        for (int i = k + 1; i < N; i++) {
            double factor = augmented[i][k] / augmented[k][k];
            for (int j = k; j <= N; j++) {
                augmented[i][j] -= factor * augmented[k][j];
            }
        }
    }
    
    U.resize(N);
    for (int i = N - 1; i >= 0; i--) {
        U[i] = augmented[i][N];
        for (int j = i + 1; j < N; j++) {
            U[i] -= augmented[i][j] * U[j];
        }
        U[i] /= augmented[i][i];
    }
}

void FDDiscretization::printSystem() const {
    std::cout << "Linear System A U = F (N = " << N << "):" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    
    for (int i = 0; i < std::min(10, N); i++) {
        for (int j = 0; j < std::min(10, N); j++) {
            std::cout << std::setw(10) << A[i][j] << " ";
        }
        if (N > 10) std::cout << "... ";
        std::cout << " | " << std::setw(10) << F[i] << std::endl;
    }
}

void FDDiscretization::computeError(double& L1, double& L2, double& Linf) const {
    L1 = 0.0;
    L2 = 0.0;
    Linf = 0.0;
    
    auto eq_points = grid.getEquationPoints();
    double h = grid.getHx();
    
    for (const auto& point : eq_points) {
        int i = point.first;
        int j = point.second;
        int idx = grid.getGlobalIdx(i, j);
        
        if (idx >= 0) {
            double x = grid.getX(i, j);
            double y = grid.getY(i, j);
            double exact_val = test_func.exact(x, y);
            double error = std::abs(U[idx] - exact_val);
            
            L1 += error * h;
            L2 += error * error * h;
            Linf = std::max(Linf, error);
        }
    }
    
    L2 = std::sqrt(L2);
}