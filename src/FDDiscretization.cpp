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

FDDiscretization::FDDiscretization(const Grid& grid_, const TestFunction& test_func_,
                                   const BoundaryConfig& left_, const BoundaryConfig& right_,
                                   const BoundaryConfig& bottom_, const BoundaryConfig& top_,
                                   const BoundaryConfig& hole_)
    : grid(grid_), test_func(test_func_),
      bc_left(left_), bc_right(right_), bc_bottom(bottom_), bc_top(top_),
      bc_hole(hole_) {
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
    
    // 正方形边界
    if (i == 0) is_dirichlet = is_dirichlet || isDirichlet("bottom");
    if (i == grid.getNy() + 1) is_dirichlet = is_dirichlet || isDirichlet("top");
    if (j == 0) is_dirichlet = is_dirichlet || isDirichlet("left");
    if (j == grid.getNx() + 1) is_dirichlet = is_dirichlet || isDirichlet("right");
    
    // 圆孔边界：如果点在圆周上且圆孔是 Dirichlet
    if (grid.hasHole() && isHoleDirichlet()) {
        double x = grid.getX(i, j);
        double y = grid.getY(i, j);
        double cx = grid.getHoleCx();
        double cy = grid.getHoleCy();
        double r = grid.getHoleR();
        double dx = x - cx;
        double dy = y - cy;
        double dist_sq = dx*dx + dy*dy;
        double eps = 1e-10;
        if (std::abs(dist_sq - r*r) < eps) {
            is_dirichlet = true;
        }
    }
    
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
            if (isDirichletOnPoint(i, j)) {
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
    double g = 0.0;
    
    // 检查是否是圆孔边界
    if (grid.hasHole() && isHoleDirichlet()) {
        double x = grid.getX(i, j);
        double y = grid.getY(i, j);
        double cx = grid.getHoleCx();
        double cy = grid.getHoleCy();
        double r = grid.getHoleR();
        double dx = x - cx;
        double dy = y - cy;
        double dist_sq = dx*dx + dy*dy;
        double eps = 1e-10;
        if (std::abs(dist_sq - r*r) < eps) {
            // 圆孔边界上的点
            g = getHoleBoundaryValue(x, y);
            A[idx][idx] = 1.0;
            F[idx] = g;
            return;
        }
    }
    
    // 正方形边界
    std::string side;
    if (i == 0) side = "bottom";
    else if (i == grid.getNy() + 1) side = "top";
    else if (j == 0) side = "left";
    else side = "right";
    
    g = getBoundaryValue(i, j, side);
    
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
        // 上下边界：内部邻居系数为 -2/h²
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
    
    // 计算法线方向（用于 Neumann 条件）
    double cx = grid.getHoleCx();
    double cy = grid.getHoleCy();
    double dx_c = x - cx;
    double dy_c = y - cy;
    double dist_c = std::sqrt(dx_c*dx_c + dy_c*dy_c);
    double nx = dx_c / dist_c;
    double ny = dy_c / dist_c;
    
    // 获取圆孔 Neumann 值（如果需要）
    double g_hole = 0.0;
    if (isHoleNeumann()) {
        g_hole = getHoleNeumannValue(x, y);
    }
    
    // ==================== Dirichlet 分支 ====================
    if (isHoleDirichlet()) {
        
        // x 方向（西边和东边）
        bool west_hole = (grid.getPointType(i, j-1) == PointType::HOLE);
        bool east_hole = (grid.getPointType(i, j+1) == PointType::HOLE);
        
        if (west_hole || east_hole) {
            double theta_w = 0.0, theta_e = 0.0;
            double u_w = 0.0, u_e = 0.0;
            
            if (west_hole) {
                double dist_w = findDistanceToHoleBoundary(i, j, "west");
                theta_w = dist_w / h;
                u_w = getHoleBoundaryValue(x - dist_w, y);
            }
            if (east_hole) {
                double dist_e = findDistanceToHoleBoundary(i, j, "east");
                theta_e = dist_e / h;
                u_e = getHoleBoundaryValue(x + dist_e, y);
            }
            
            if (west_hole) {
                // 西边有 hole
                double coeff_x = 2.0 / (theta_w * (1.0 + theta_w) * h2);
                center_coeff += (1.0 + theta_w) * coeff_x;
                rhs += coeff_x * u_w;
                // 东边是正常点，作为未知量
                int east_idx = grid.getGlobalIdx(i, j+1);
                if (east_idx >= 0) {
                    A[idx][east_idx] = -theta_w * coeff_x;
                }
            } else if (east_hole) {
                // 东边有 hole
                double coeff_x = 2.0 / (theta_e * (1.0 + theta_e) * h2);
                center_coeff += (1.0 + theta_e) * coeff_x;
                rhs += coeff_x * u_e;
                // 西边是正常点，作为未知量
                int west_idx = grid.getGlobalIdx(i, j-1);
                if (west_idx >= 0) {
                    A[idx][west_idx] = -theta_e * coeff_x;
                }
            }
        } else {
            // x 方向没有 hole，正常处理
            int west_idx = grid.getGlobalIdx(i, j-1);
            if (west_idx >= 0) {
                A[idx][west_idx] = -1.0 / h2;
                center_coeff += 1.0 / h2;
            } else if (grid.getPointType(i, j-1) == PointType::BOUNDARY && isDirichletOnPoint(i, j-1)) {
                double bc_val = getBoundaryValue(i, j-1, "left");
                rhs += bc_val / h2;
                center_coeff += 1.0 / h2;
            }
            
            int east_idx = grid.getGlobalIdx(i, j+1);
            if (east_idx >= 0) {
                A[idx][east_idx] = -1.0 / h2;
                center_coeff += 1.0 / h2;
            } else if (grid.getPointType(i, j+1) == PointType::BOUNDARY && isDirichletOnPoint(i, j+1)) {
                double bc_val = getBoundaryValue(i, j+1, "right");
                rhs += bc_val / h2;
                center_coeff += 1.0 / h2;
            }
        }
        
        // y 方向（南边和北边）
        bool south_hole = (grid.getPointType(i-1, j) == PointType::HOLE);
        bool north_hole = (grid.getPointType(i+1, j) == PointType::HOLE);
        
        if (south_hole || north_hole) {
            double alpha_s = 0.0, alpha_n = 0.0;
            double u_s = 0.0, u_n = 0.0;
            
            if (south_hole) {
                double dist_s = findDistanceToHoleBoundary(i, j, "south");
                alpha_s = dist_s / h;
                u_s = getHoleBoundaryValue(x, y - dist_s);
            }
            if (north_hole) {
                double dist_n = findDistanceToHoleBoundary(i, j, "north");
                alpha_n = dist_n / h;
                u_n = getHoleBoundaryValue(x, y + dist_n);
            }
            
            if (south_hole) {
                // 南边有 hole
                double coeff_y = 2.0 / (alpha_s * (1.0 + alpha_s) * h2);
                center_coeff += (1.0 + alpha_s) * coeff_y;
                rhs += coeff_y * u_s;
                // 北边是正常点，作为未知量
                int north_idx = grid.getGlobalIdx(i+1, j);
                if (north_idx >= 0) {
                    A[idx][north_idx] = -alpha_s * coeff_y;
                }
            } else if (north_hole) {
                // 北边有 hole
                double coeff_y = 2.0 / (alpha_n * (1.0 + alpha_n) * h2);
                center_coeff += (1.0 + alpha_n) * coeff_y;
                rhs += coeff_y * u_n;
                // 南边是正常点
                int south_idx = grid.getGlobalIdx(i-1, j);
                if (south_idx >= 0) {
                    A[idx][south_idx] = -alpha_n * coeff_y;
                }
            }
        } else {
            // y 方向没有 hole，正常处理
            int south_idx = grid.getGlobalIdx(i-1, j);
            if (south_idx >= 0) {
                A[idx][south_idx] = -1.0 / h2;
                center_coeff += 1.0 / h2;
            } else if (grid.getPointType(i-1, j) == PointType::BOUNDARY && isDirichletOnPoint(i-1, j)) {
                double bc_val = getBoundaryValue(i-1, j, "bottom");
                rhs += bc_val / h2;
                center_coeff += 1.0 / h2;
            }
            
            int north_idx = grid.getGlobalIdx(i+1, j);
            if (north_idx >= 0) {
                A[idx][north_idx] = -1.0 / h2;
                center_coeff += 1.0 / h2;
            } else if (grid.getPointType(i+1, j) == PointType::BOUNDARY && isDirichletOnPoint(i+1, j)) {
                double bc_val = getBoundaryValue(i+1, j, "top");
                rhs += bc_val / h2;
                center_coeff += 1.0 / h2;
            }
        }
    }
    
    // ==================== Neumann 分支 ====================
    else if (isHoleNeumann()) {
        
        // 检查哪些方向有 HOLE
        bool west_hole = (grid.getPointType(i, j-1) == PointType::HOLE);
        bool east_hole = (grid.getPointType(i, j+1) == PointType::HOLE);
        bool south_hole = (grid.getPointType(i-1, j) == PointType::HOLE);
        bool north_hole = (grid.getPointType(i+1, j) == PointType::HOLE);
        int south_idx = grid.getGlobalIdx(i-1, j);
        int north_idx = grid.getGlobalIdx(i+1, j);
        int west_idx = grid.getGlobalIdx(i, j-1);
        int east_idx = grid.getGlobalIdx(i, j+1);
        rhs = h * g_hole;
        
        if (north_hole && east_hole) {
            A[idx][south_idx] = -ny;
            A[idx][west_idx] = -nx;
            center_coeff = nx + ny;
            }
        else if (south_hole && west_hole) {
            A[idx][north_idx] = ny;
            A[idx][east_idx] = nx;
            center_coeff = -nx + -ny;
        }
        else if (west_hole && north_hole) {
            A[idx][east_idx] = nx;
            A[idx][south_idx] = -ny;
            center_coeff = -nx + ny;
        }    
        else if (east_hole && south_hole){
            A[idx][west_idx] = -nx;
            A[idx][north_idx] = ny;
            center_coeff = nx + -ny;
        }
        else if (east_hole){
            A[idx][west_idx] = -nx;
            A[idx][north_idx] = ny;
            center_coeff = nx + -ny;            
        }
        else if (south_hole) {
            A[idx][north_idx] = ny;
            A[idx][east_idx] = nx;
            center_coeff = -nx + -ny;
        }
        else if (west_hole) {
            A[idx][east_idx] = nx;
            A[idx][south_idx] = -ny;
            center_coeff = -nx + ny;
        }    
        else{
            A[idx][south_idx] = -ny;
            A[idx][west_idx] = -nx;
            center_coeff = nx + ny;;
        }

    }
    
    A[idx][idx] = center_coeff;
    F[idx] = rhs;
}

double FDDiscretization::findDistanceToHoleBoundary(int i, int j, const std::string& direction) {
    double x = grid.getX(i, j);
    double y = grid.getY(i, j);
    double h = grid.getHx();
    
    // 从 Grid 获取圆孔参数
    double cx = grid.getHoleCx();
    double cy = grid.getHoleCy();
    double r = grid.getHoleR();
    
    // 解析计算沿给定方向到圆的距离
    double dx_dir = 0.0, dy_dir = 0.0;
    if (direction == "west") { dx_dir = -1.0; dy_dir = 0.0; }
    else if (direction == "east") { dx_dir = 1.0; dy_dir = 0.0; }
    else if (direction == "south") { dx_dir = 0.0; dy_dir = -1.0; }
    else if (direction == "north") { dx_dir = 0.0; dy_dir = 1.0; }
    
    double dx = x - cx;
    double dy = y - cy;
    
    // 二次方程系数: a*t^2 + b*t + c = 0
    double a = dx_dir*dx_dir + dy_dir*dy_dir;  // = 1
    double b = 2 * (dx*dx_dir + dy*dy_dir);
    double c = dx*dx + dy*dy - r*r;
    
    double discriminant = b*b - 4*a*c;
    if (discriminant < 0) {
        return h;
    }
    
    double sqrt_disc = std::sqrt(discriminant);
    double t1 = (-b - sqrt_disc) / (2*a);
    double t2 = (-b + sqrt_disc) / (2*a);
    
    // 选择正的小的 t
    double t = h;
    if (t1 > 1e-12 && t1 < t) t = t1;
    if (t2 > 1e-12 && t2 < t) t = t2;
    
    return t;
}

double FDDiscretization::getHoleNeumannValue(double x, double y) const {
    // 从 Grid 获取圆孔参数
    double cx = grid.getHoleCx();
    double cy = grid.getHoleCy();
    double r = grid.getHoleR();
    
    // 计算从圆心到 P 的方向（法线方向）
    double dx = x - cx;
    double dy = y - cy;
    double dist = std::sqrt(dx*dx + dy*dy);
    
    // 单位法向量（从圆心指向外）
    double nx = dx / dist;
    double ny = dy / dist;
    
    // 调用 TestFunction 的圆孔 Neumann 计算函数
    return test_func.getHoleNeumannBC(x, y, nx, ny);
}



// 获取圆孔边界值（Dirichlet）
double FDDiscretization::getHoleBoundaryValue(double x, double y) const {
    return test_func.exact(x, y);
}

// 判断圆孔边界是否为 Dirichlet
bool FDDiscretization::isHoleDirichlet() const {
    return bc_hole.type == "dirichlet";
}

// 判断圆孔边界是否为 Neumann
bool FDDiscretization::isHoleNeumann() const {
    return bc_hole.type == "neumann";
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
    
    // 分类统计变量
    double L1_regular = 0.0, L1_irregular = 0.0, L1_boundary = 0.0;
    double L2_regular = 0.0, L2_irregular = 0.0, L2_boundary = 0.0;
    double Linf_regular = 0.0, Linf_irregular = 0.0, Linf_boundary = 0.0;
    int count_regular = 0, count_irregular = 0, count_boundary = 0;
    
    auto eq_points = grid.getEquationPoints();
    double hx = grid.getHx();
    double hy = grid.getHy();
    double area = hx * hy;
    
    for (const auto& point : eq_points) {
        int i = point.first;
        int j = point.second;
        int idx = grid.getGlobalIdx(i, j);
        
        if (idx >= 0 && idx < (int)U.size()) {
            double x = grid.getX(i, j);
            double y = grid.getY(i, j);
            double exact_val = test_func.exact(x, y);
            double error = std::abs(U[idx] - exact_val);
            
            // 总体误差
            L1 += error * area;
            L2 += error * error * area;
            Linf = std::max(Linf, error);
            
            // 根据点类型分类统计
            PointType type = grid.getPointType(i, j);
            if (type == PointType::REGULAR) {
                L1_regular += error * area;
                L2_regular += error * error * area;
                Linf_regular = std::max(Linf_regular, error);
                count_regular++;
            } else if (type == PointType::IRREGULAR) {
                L1_irregular += error * area;
                L2_irregular += error * error * area;
                Linf_irregular = std::max(Linf_irregular, error);
                count_irregular++;
            } else if (type == PointType::BOUNDARY) {
                L1_boundary += error * area;
                L2_boundary += error * error * area;
                Linf_boundary = std::max(Linf_boundary, error);
                count_boundary++;
            }
        }
    }
    
    L2 = std::sqrt(L2);
    
    // 计算分类的 L2
    L2_regular = std::sqrt(L2_regular);
    L2_irregular = std::sqrt(L2_irregular);
    L2_boundary = std::sqrt(L2_boundary);
    
    // 打印分类统计结果
    std::cout << "\n=== Detailed Error Analysis ===" << std::endl;
    std::cout << "Total points: " << eq_points.size() << std::endl;
    std::cout << "  Regular points:   " << count_regular 
              << " (L1=" << L1_regular << ", L2=" << L2_regular << ", Linf=" << Linf_regular << ")" << std::endl;
    std::cout << "  Irregular points: " << count_irregular 
              << " (L1=" << L1_irregular << ", L2=" << L2_irregular << ", Linf=" << Linf_irregular << ")" << std::endl;
    std::cout << "  Boundary points:  " << count_boundary 
              << " (L1=" << L1_boundary << ", L2=" << L2_boundary << ", Linf=" << Linf_boundary << ")" << std::endl;
    
    // 计算平均误差（每个点）
    if (count_regular > 0) {
        std::cout << "  Average regular error: " << (L1_regular / area / count_regular) << std::endl;
    }
    if (count_irregular > 0) {
        std::cout << "  Average irregular error: " << (L1_irregular / area / count_irregular) << std::endl;
    }
    if (count_boundary > 0) {
        std::cout << "  Average boundary error: " << (L1_boundary / area / count_boundary) << std::endl;
    }
}