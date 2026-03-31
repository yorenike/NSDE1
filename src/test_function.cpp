#include <iostream>
#include <iomanip>
#include "TestFunction.h"

int main() {
    std::vector<std::string> test_names = {"exp_y_sin_x", "polynomial", "trigonometric"};
    
    for (const auto& name : test_names) {
        std::cout << "========== Testing: " << name << " ==========" << std::endl;
        
        try {
            TestFunction tf(name);
            
            // 测试点
            double x = 0.25;
            double y = 0.5;
            
            std::cout << std::fixed << std::setprecision(6);
            std::cout << "At (x=" << x << ", y=" << y << "):" << std::endl;
            std::cout << "  u(x,y) = " << tf.exact(x, y) << std::endl;
            std::cout << "  f(x,y) = " << tf.f(x, y) << std::endl;
            
            std::cout << "\nBoundary conditions at (x=" << x << ", y=" << y << "):" << std::endl;
            std::cout << "  Left (x=0):   " << tf.getDirichletBC(0, y, "left") << std::endl;
            std::cout << "  Right (x=1):  " << tf.getDirichletBC(1, y, "right") << std::endl;
            std::cout << "  Bottom (y=0): " << tf.getDirichletBC(x, 0, "bottom") << std::endl;
            std::cout << "  Top (y=1):    " << tf.getDirichletBC(x, 1, "top") << std::endl;
            
            std::cout << "\nNeumann BC (∂u/∂n):" << std::endl;
            std::cout << "  Left:   " << tf.getNeumannBC(0, y, "left") << std::endl;
            std::cout << "  Right:  " << tf.getNeumannBC(1, y, "right") << std::endl;
            std::cout << "  Bottom: " << tf.getNeumannBC(x, 0, "bottom") << std::endl;
            std::cout << "  Top:    " << tf.getNeumannBC(x, 1, "top") << std::endl;
            
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
        
        std::cout << std::endl;
    }
    
    return 0;
}