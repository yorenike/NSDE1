#ifndef TEST_FUNCTION_H
#define TEST_FUNCTION_H

#include <string>
#include <functional>
#include <cmath>
#include <stdexcept>

class TestFunction {
private:
    std::string name;
    
    // 解析解 u(x,y)
    std::function<double(double, double)> u_func;
    
    // 右端项 f(x,y) = -Δu
    std::function<double(double, double)> f_func;
    
    // Dirichlet 边界条件
    std::function<double(double, double)> bc_left;
    std::function<double(double, double)> bc_right;
    std::function<double(double, double)> bc_bottom;
    std::function<double(double, double)> bc_top;
    
    // Neumann 边界条件（法向导数）
    std::function<double(double, double)> du_dn_left;
    std::function<double(double, double)> du_dn_right;
    std::function<double(double, double)> du_dn_bottom;
    std::function<double(double, double)> du_dn_top;
    
    // 初始化各个测试函数
    void initExpYSinX();
    void initPolynomial();
    void initTrigonometric();
    
public:
    // 构造函数：根据名称选择测试函数
    TestFunction(const std::string& name);
    
    // 获取解析解
    double exact(double x, double y) const;
    
    // 获取右端项 f = -Δu
    double f(double x, double y) const;
    
    // 获取 Dirichlet 边界条件值
    double getDirichletBC(double x, double y, const std::string& side) const;
    
    // 获取 Neumann 边界条件值（法向导数）
    double getNeumannBC(double x, double y, const std::string& side) const;
    
    // 获取测试函数名称
    std::string getName() const { return name; }
    
    // 判断边界值是否来自测试函数
    static bool isFromTestFunction(const std::string& value) {
        return value == "from_test_function";
    }
};

#endif // TEST_FUNCTION_H