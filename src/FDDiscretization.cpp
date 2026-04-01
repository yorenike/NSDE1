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

bool FDDiscretization::isNeumann(const std::string& side) const {
    if (side == "left") return bc_left.type == "neumann";
    if (side == "right") return bc_right.type == "neumann";
    if (side == "bottom") return bc_bottom.type == "neumann";
    if (side == "top") return bc_top.type == "neumann";
    return false;
}

bool FDDiscretization::isDirichlet(const std::string& side) const {
    if (side == "left") return bc_left.type == "dirichlet";
    if (side == "right") return bc_right.type == "dirichlet";
    if (side == "bottom") return bc_bottom.type == "dirichlet";
    if (side == "top") return bc_top.type == "dirichlet";
    return false;
}

bool FDDiscretization::isDirichletOnPoint(int i, int j) const {
    bool is_dirichlet = false;
    
    if (i == 0) is_dirichlet = is_dirichlet || isDirichlet("bottom");
    if (i == grid.getNy() + 1) is_dirichlet = is_dirichlet || isDirichlet("top");
    if (j == 0) is_dirichlet = is_dirichlet || isDirichlet("left");
    if (j == grid.getNx() + 1) is_dirichlet = is_dirichlet || isDirichlet("right");
    
    return is_dirichlet;
}

bool FDDiscretization::isNeumannOnPoint(int i, int j) const {
    bool is_neumann = false;
    
    if (i == 0) is_neumann = is_neumann || isNeumann("bottom");
    if (i == grid.getNy() + 1) is_neumann = is_neumann || isNeumann("top");
    if (j == 0) is_neumann = is_neumann || isNeumann("left");
    if (j == grid.getNx() + 1) is_neumann = is_neumann || isNeumann("right");
    
    return is_neumann;
}

bool FDDiscretization::isCornerPoint(int i, int j) const {
    return (i == 0 || i == grid.getNy() + 1) && (j == 0 || j == grid.getNx() + 1);
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
        return 0.0;
    }
}

double FDDiscretization::getNeumannValue(int i, int j, const std::string& side) const {
    double x = grid.getX(i, j);
    double y = grid.getY(i, j);
    
    const BoundaryConfig* bc = nullptr;
    if (side == "left") bc = &bc_left;
    else if (side == "right") bc = &bc_right;
    else if (side == "bottom") bc = &bc_bottom;
    else if (side == "top") bc = &bc_top;
    else return 0.0;
    
    if (bc->from_test) {
        return test_func.getNeumannBC(x, y, side);
    } else {
        return 0.0;
    }
}

void FDDiscretization::assemble() {
    initializeSystem();
    
    double h = grid.getHx();
    double h2 = h * h;
    
    for (int i = 0; i <= grid.getNy() + 1; i++) {
        for (int j = 0; j <= grid.getNx() + 1; j++) {
            int idx = grid.getGlobalIdx(i, j);
            if (idx < 0) continue;
            
            PointType type = grid.getPointType(i, j);
            
            // 优先处理 Dirichlet：只要任意一个边界是 Dirichlet，该点就是 Dirichlet
            if (type == PointType::BOUNDARY && isDirichletOnPoint(i, j)) {
                assembleDirichletPoint(i, j, idx);
                continue;
            }
            
            // 处理两个 Neumann 的角点
            if (type == PointType::BOUNDARY && isCornerPoint(i, j) && isNeumannOnPoint(i, j)) {
                assembleNeumannCornerPoint(i, j, idx);
                continue;
            }
            
            // 处理 Neumann 边界点（非角点）
            if (type == PointType::BOUNDARY && isNeumannOnPoint(i, j)) {
                assembleNeumannBoundaryPoint(i, j, idx);
                continue;
            }
            
            // 处理内部点
            if (type == PointType::REGULAR) {
                assembleRegularPoint(i, j, idx);
            } else if (type == PointType::IRREGULAR) {
                assembleIrregularPoint(i, j, idx);
            }
        }
    }
}

void FDDiscretization::assembleDirichletPoint(int i, int j, int idx) {
    std::string side;
    if (i == 0) side = "bottom";
    else if (i == grid.getNy() + 1) side = "top";
    else if (j == 0) side = "left";
    else side = "right";
    
    double g = getBoundaryValue(i, j, side);
    
    A[idx][idx] = 1.0;
    F[idx] = g;
}

void FDDiscretization::assembleNeumannBoundaryPoint(int i, int j, int idx) {
    double h = grid.getHx();
    double h2 = h * h;
    double x = grid.getX(i, j);
    double y = grid.getY(i, j);
    
    std::string side;
    if (i == 0) side = "bottom";
    else if (i == grid.getNy() + 1) side = "top";
    else if (j == 0) side = "left";
    else side = "right";
    
    double sigma = getNeumannValue(i, j, side);
    
    // 设置中心点系数
    A[idx][idx] = 4.0 / h2;
    
    // 添加内部邻居的贡献（区分边界类型）
    if (side == "left" || side == "right") {
        // 左右边界：内部邻居系数为 -2/h²
        int i_in = i;
        int j_in = (side == "left") ? j + 1 : j - 1;
        int idx_in = grid.getGlobalIdx(i_in, j_in);
        if (idx_in >= 0) {
            A[idx][idx_in] -= 2.0 / h2;
        }
    } else {
        // 上下边界：内部邻居系数为 -1/h²
        int i_in = (side == "bottom") ? i + 1 : i - 1;
        int j_in = j;
        int idx_in = grid.getGlobalIdx(i_in, j_in);
        if (idx_in >= 0) {
            A[idx][idx_in] -= 2.0 / h2;
        }
    }
    
    // 添加切线方向的邻居（系数 -1/h²）
    if (side == "left" || side == "right") {
        if (i > 0) {
            int idx_down = grid.getGlobalIdx(i-1, j);
            if (idx_down >= 0) A[idx][idx_down] = -1.0 / h2;
        }
        if (i < grid.getNy() + 1) {
            int idx_up = grid.getGlobalIdx(i+1, j);
            if (idx_up >= 0) A[idx][idx_up] = -1.0 / h2;
        }
    } else {
        if (j > 0) {
            int idx_left = grid.getGlobalIdx(i, j-1);
            if (idx_left >= 0) A[idx][idx_left] = -1.0 / h2;
        }
        if (j < grid.getNx() + 1) {
            int idx_right = grid.getGlobalIdx(i, j+1);
            if (idx_right >= 0) A[idx][idx_right] = -1.0 / h2;
        }
    }
    
    F[idx] = test_func.f(x, y) + 2.0 * sigma / h;
}

void FDDiscretization::assembleNeumannCornerPoint(int i, int j, int idx) {
    double h = grid.getHx();
    double h2 = h * h;
    double x = grid.getX(i, j);
    double y = grid.getY(i, j);
    
    // 确定两个边界的方向和 Neumann 值
    double sigma1 = 0.0, sigma2 = 0.0;
    
    if (i == 0) {
        sigma1 = getNeumannValue(i, j, "bottom");
    } else if (i == grid.getNy() + 1) {
        sigma1 = getNeumannValue(i, j, "top");
    }
    
    if (j == 0) {
        sigma2 = getNeumannValue(i, j, "left");
    } else if (j == grid.getNx() + 1) {
        sigma2 = getNeumannValue(i, j, "right");
    }
    
    // 确定内部邻居
    int i_in = (i == 0) ? i + 1 : i - 1;
    int j_in = (j == 0) ? j + 1 : j - 1;
    
    int idx_in_i = grid.getGlobalIdx(i_in, j);
    int idx_in_j = grid.getGlobalIdx(i, j_in);
    
    // 角点 Neumann 公式
    // -2U_{i_in,j} - 2U_{i,j_in} + 4U_{i,j} = h²f + 2h(σ1 + σ2)
    A[idx][idx] = 4.0 / h2;
    
    if (idx_in_i >= 0) {
        A[idx][idx_in_i] = -2.0 / h2;
    }
    if (idx_in_j >= 0) {
        A[idx][idx_in_j] = -2.0 / h2;
    }
    
    // 右端项
    F[idx] = test_func.f(x, y) + 2.0 * (sigma1 + sigma2) / h;
}

void FDDiscretization::assembleRegularPoint(int i, int j, int idx) {
    double h = grid.getHx();
    double h2 = h * h;
    double x = grid.getX(i, j);
    double y = grid.getY(i, j);
    
    A[idx][idx] = 4.0 / h2;
    
    // 西边
    int neighbor_idx = grid.getGlobalIdx(i, j-1);
    if (neighbor_idx >= 0) {
        A[idx][neighbor_idx] = -1.0 / h2;
    } else if (grid.getPointType(i, j-1) == PointType::BOUNDARY && isDirichletOnPoint(i, j-1)) {
        double bc_val = getBoundaryValue(i, j-1, "left");
        F[idx] += bc_val / h2;
    }
    
    // 东边
    neighbor_idx = grid.getGlobalIdx(i, j+1);
    if (neighbor_idx >= 0) {
        A[idx][neighbor_idx] = -1.0 / h2;
    } else if (grid.getPointType(i, j+1) == PointType::BOUNDARY && isDirichletOnPoint(i, j+1)) {
        double bc_val = getBoundaryValue(i, j+1, "right");
        F[idx] += bc_val / h2;
    }
    
    // 南边
    neighbor_idx = grid.getGlobalIdx(i-1, j);
    if (neighbor_idx >= 0) {
        A[idx][neighbor_idx] = -1.0 / h2;
    } else if (grid.getPointType(i-1, j) == PointType::BOUNDARY && isDirichletOnPoint(i-1, j)) {
        double bc_val = getBoundaryValue(i-1, j, "bottom");
        F[idx] += bc_val / h2;
    }
    
    // 北边
    neighbor_idx = grid.getGlobalIdx(i+1, j);
    if (neighbor_idx >= 0) {
        A[idx][neighbor_idx] = -1.0 / h2;
    } else if (grid.getPointType(i+1, j) == PointType::BOUNDARY && isDirichletOnPoint(i+1, j)) {
        double bc_val = getBoundaryValue(i+1, j, "top");
        F[idx] += bc_val / h2;
    }
    
    // 右端项：f(x,y)
    F[idx] += test_func.f(x, y);
}

void FDDiscretization::assembleIrregularPoint(int i, int j, int idx) {
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
        double x_boundary = x - dist;
        double y_boundary = y;
        rhs += coeff * theta * test_func.exact(x_boundary, y_boundary);
    } else {
        int neighbor_idx = grid.getGlobalIdx(i, j-1);
        if (neighbor_idx >= 0) {
            A[idx][neighbor_idx] = -1.0 / h2;
            center_coeff += 1.0 / h2;
        } else if (grid.getPointType(i, j-1) == PointType::BOUNDARY && isDirichletOnPoint(i, j-1)) {
            double bc_val = getBoundaryValue(i, j-1, "left");
            rhs += bc_val / h2;
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
    } else {
        int neighbor_idx = grid.getGlobalIdx(i, j+1);
        if (neighbor_idx >= 0) {
            A[idx][neighbor_idx] = -1.0 / h2;
            center_coeff += 1.0 / h2;
        } else if (grid.getPointType(i, j+1) == PointType::BOUNDARY && isDirichletOnPoint(i, j+1)) {
            double bc_val = getBoundaryValue(i, j+1, "right");
            rhs += bc_val / h2;
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
    } else {
        int neighbor_idx = grid.getGlobalIdx(i-1, j);
        if (neighbor_idx >= 0) {
            A[idx][neighbor_idx] = -1.0 / h2;
            center_coeff += 1.0 / h2;
        } else if (grid.getPointType(i-1, j) == PointType::BOUNDARY && isDirichletOnPoint(i-1, j)) {
            double bc_val = getBoundaryValue(i-1, j, "bottom");
            rhs += bc_val / h2;
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
    } else {
        int neighbor_idx = grid.getGlobalIdx(i+1, j);
        if (neighbor_idx >= 0) {
            A[idx][neighbor_idx] = -1.0 / h2;
            center_coeff += 1.0 / h2;
        } else if (grid.getPointType(i+1, j) == PointType::BOUNDARY && isDirichletOnPoint(i+1, j)) {
            double bc_val = getBoundaryValue(i+1, j, "top");
            rhs += bc_val / h2;
            center_coeff += 1.0 / h2;
        }
    }
    
    A[idx][idx] = center_coeff;
    F[idx] = rhs;
}

double FDDiscretization::findDistanceToHoleBoundary(int i, int j, const std::string& direction) {
    double x = grid.getX(i, j);
    double y = grid.getY(i, j);
    double h = grid.getHx();
    
    double cx = grid.getHoleCx();
    double cy = grid.getHoleCy();
    double r = grid.getHoleR();
    
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
    double hx = grid.getHx();
    double hy = grid.getHy();
    
    for (const auto& point : eq_points) {
        int i = point.first;
        int j = point.second;
        int idx = grid.getGlobalIdx(i, j);
        
        if (idx >= 0 && idx < (int)U.size()) {
            double x = grid.getX(i, j);
            double y = grid.getY(i, j);
            double exact_val = test_func.exact(x, y);
            double error = std::abs(U[idx] - exact_val);
            
            L1 += error * hx * hy;
            L2 += error * error * hx * hy;
            Linf = std::max(Linf, error);
        }
    }
    
    L2 = std::sqrt(L2);
}