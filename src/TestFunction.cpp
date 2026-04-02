#include "TestFunction.h"

// ========== 测试函数1: exp(y + sin(x)) ==========
// 这是题目要求的测试函数
void TestFunction::initExpYSinX() {
    name = "exp_y_sin_x";
    
    // 解析解
    u_func = [](double x, double y) {
        return std::exp(y + std::sin(x));
    };
    
    // 右端项 f = -Δu
    // u = exp(y + sin(x))
    // u_x = exp(y + sin(x)) * cos(x)
    // u_xx = exp(y + sin(x)) * (cos^2(x) - sin(x))
    // u_y = exp(y + sin(x))
    // u_yy = exp(y + sin(x))
    // Δu = exp(y + sin(x)) * (cos^2(x) - sin(x) + 1)
    // f = -Δu = -exp(y + sin(x)) * (cos^2(x) - sin(x) + 1)
    f_func = [](double x, double y) {
        double u = std::exp(y + std::sin(x));
        double cosx = std::cos(x);
        return -u * (cosx * cosx - std::sin(x) + 1.0);
    };
    
    // Dirichlet 边界条件
    bc_left = [](double x, double y) { return std::exp(y + std::sin(0)); };      // x=0
    bc_right = [](double x, double y) { return std::exp(y + std::sin(1)); };     // x=1
    bc_bottom = [](double x, double y) { return std::exp(0 + std::sin(x)); };    // y=0
    bc_top = [](double x, double y) { return std::exp(1 + std::sin(x)); };       // y=1
    
    // Neumann 边界条件（法向导数）
    // ∂u/∂n = ∇u · n
    // 左边界 (x=0, n=(-1,0)): ∂u/∂n = -u_x = -u * cos(0) = -u
    // 右边界 (x=1, n=(1,0)): ∂u/∂n = u_x = u * cos(1)
    // 下边界 (y=0, n=(0,-1)): ∂u/∂n = -u_y = -u
    // 上边界 (y=1, n=(0,1)): ∂u/∂n = u_y = u
    du_dn_left = [](double x, double y) {
        double u = std::exp(y + std::sin(x));
        return -u * std::cos(x);
    };
    du_dn_right = [](double x, double y) {
        double u = std::exp(y + std::sin(x));
        return u * std::cos(x);
    };
    du_dn_bottom = [](double x, double y) {
        double u = std::exp(y + std::sin(x));
        return -u;
    };
    du_dn_top = [](double x, double y) {
        double u = std::exp(y + std::sin(x));
        return u;
    };
}

// ========== 测试函数2: 多项式 u(x,y) = x(1-x)y(1-y) ==========
// 这是一个简单的多项式函数，满足齐次 Dirichlet 边界条件
void TestFunction::initPolynomial() {
    name = "polynomial";
    
    // 解析解
    u_func = [](double x, double y) {
        return x * (1.0 - x) * y * (1.0 - y);
    };
    
    // 右端项 f = -Δu
    // u = (x - x^2)(y - y^2)
    // u_xx = -2 * (y - y^2)
    // u_yy = -2 * (x - x^2)
    // Δu = -2(y - y^2) - 2(x - x^2) = -2(x + y - x^2 - y^2)
    // f = -Δu = 2(x + y - x^2 - y^2)
    f_func = [](double x, double y) {
        return 2.0 * (x + y - x*x - y*y);
    };
    
    // Dirichlet 边界条件：在边界上 u=0
    bc_left = [](double x, double y) { return 0.0; };
    bc_right = [](double x, double y) { return 0.0; };
    bc_bottom = [](double x, double y) { return 0.0; };
    bc_top = [](double x, double y) { return 0.0; };
    
    // Neumann 边界条件（法向导数）
    // u = x(1-x)y(1-y)
    // u_x = (1-2x)y(1-y)
    // u_y = x(1-x)(1-2y)
    du_dn_left = [](double x, double y) {
        // 左边界 x=0, n=(-1,0): ∂u/∂n = -u_x = -[(1-2*0)y(1-y)] = -y(1-y)
        return -y * (1.0 - y);
    };
    du_dn_right = [](double x, double y) {
        // 右边界 x=1, n=(1,0): ∂u/∂n = u_x = (1-2*1)y(1-y) = -y(1-y)
        return -y * (1.0 - y);
    };
    du_dn_bottom = [](double x, double y) {
        // 下边界 y=0, n=(0,-1): ∂u/∂n = -u_y = -[x(1-x)(1-2*0)] = -x(1-x)
        return -x * (1.0 - x);
    };
    du_dn_top = [](double x, double y) {
        // 上边界 y=1, n=(0,1): ∂u/∂n = u_y = x(1-x)(1-2*1) = -x(1-x)
        return -x * (1.0 - x);
    };
}

// ========== 测试函数3: 三角函数 u(x,y) = sin(πx) sin(πy) ==========
// 这是一个光滑的三角函数，满足齐次 Dirichlet 边界条件
void TestFunction::initTrigonometric() {
    name = "trigonometric";
    
    // 解析解
    u_func = [](double x, double y) {
        return std::sin(M_PI * x) * std::sin(M_PI * y);
    };
    
    // 右端项 f = -Δu
    // u_xx = -π² sin(πx) sin(πy)
    // u_yy = -π² sin(πx) sin(πy)
    // Δu = -2π² sin(πx) sin(πy)
    // f = -Δu = 2π² sin(πx) sin(πy)
    f_func = [](double x, double y) {
        return 2.0 * M_PI * M_PI * std::sin(M_PI * x) * std::sin(M_PI * y);
    };
    
    // Dirichlet 边界条件：在边界上 u=0
    bc_left = [](double x, double y) { return 0.0; };
    bc_right = [](double x, double y) { return 0.0; };
    bc_bottom = [](double x, double y) { return 0.0; };
    bc_top = [](double x, double y) { return 0.0; };
    
    // Neumann 边界条件（法向导数）
    // u_x = π cos(πx) sin(πy)
    // u_y = π sin(πx) cos(πy)
    du_dn_left = [](double x, double y) {
        // 左边界 x=0, n=(-1,0): ∂u/∂n = -u_x = -π cos(0) sin(πy) = -π sin(πy)
        return -M_PI * std::sin(M_PI * y);
    };
    du_dn_right = [](double x, double y) {
        // 右边界 x=1, n=(1,0): ∂u/∂n = u_x = π cos(π) sin(πy) = -π sin(πy)
        return -M_PI * std::sin(M_PI * y);
    };
    du_dn_bottom = [](double x, double y) {
        // 下边界 y=0, n=(0,-1): ∂u/∂n = -u_y = -π sin(πx) cos(0) = -π sin(πx)
        return -M_PI * std::sin(M_PI * x);
    };
    du_dn_top = [](double x, double y) {
        // 上边界 y=1, n=(0,1): ∂u/∂n = u_y = π sin(πx) cos(π) = -π sin(πx)
        return -M_PI * std::sin(M_PI * x);
    };
}

// ========== 构造函数 ==========
TestFunction::TestFunction(const std::string& name_) {
    if (name_ == "exp_y_sin_x") {
        initExpYSinX();
    }
    else if (name_ == "polynomial") {
        initPolynomial();
    }
    else if (name_ == "trigonometric") {
        initTrigonometric();
    }
    else {
        throw std::runtime_error("Unknown test function: " + name_);
    }
}

// ========== 公共接口实现 ==========
double TestFunction::exact(double x, double y) const {
    return u_func(x, y);
}

double TestFunction::f(double x, double y) const {
    return f_func(x, y);
}

double TestFunction::getDirichletBC(double x, double y, const std::string& side) const {
    if (side == "left") return bc_left(x, y);
    if (side == "right") return bc_right(x, y);
    if (side == "bottom") return bc_bottom(x, y);
    if (side == "top") return bc_top(x, y);
    throw std::runtime_error("Unknown side: " + side);
}

double TestFunction::getNeumannBC(double x, double y, const std::string& side) const {
    if (side == "left") return du_dn_left(x, y);
    if (side == "right") return du_dn_right(x, y);
    if (side == "bottom") return du_dn_bottom(x, y);
    if (side == "top") return du_dn_top(x, y);
    throw std::runtime_error("Unknown side: " + side);
}

// ========== 获取圆孔边界上的 Neumann 条件 ==========
double TestFunction::getHoleNeumannBC(double x, double y, double nx, double ny) const {
    // 根据不同的测试函数计算法向导数
    if (name == "exp_y_sin_x") {
        // u = exp(y + sin(x))
        // ∇u = (u cos x, u)
        double u = exact(x, y);
        double cosx = std::cos(x);
        return u * (cosx * nx + ny);
    }
    else if (name == "polynomial") {
        // u = x(1-x)y(1-y)
        // u_x = (1-2x)y(1-y)
        // u_y = x(1-x)(1-2y)
        double ux = (1.0 - 2.0 * x) * y * (1.0 - y);
        double uy = x * (1.0 - x) * (1.0 - 2.0 * y);
        return ux * nx + uy * ny;
    }
    else if (name == "trigonometric") {
        // u = sin(πx) sin(πy)
        // u_x = π cos(πx) sin(πy)
        // u_y = π sin(πx) cos(πy)
        double ux = M_PI * std::cos(M_PI * x) * std::sin(M_PI * y);
        double uy = M_PI * std::sin(M_PI * x) * std::cos(M_PI * y);
        return ux * nx + uy * ny;
    }
    
    throw std::runtime_error("Unknown test function for hole Neumann BC: " + name);
}