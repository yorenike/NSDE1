#ifndef BOUNDARY_CONFIG_H
#define BOUNDARY_CONFIG_H

#include <string>

struct BoundaryConfig {
    std::string type;      // "dirichlet" 或 "neumann"
    std::string value;     // 表达式或 "from_test_function"
    bool from_test;        // 是否从测试函数获取
};

#endif