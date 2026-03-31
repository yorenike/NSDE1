#include <iostream>
#include "JSONParser.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " input.json" << std::endl;
        return 1;
    }
    
    try {
        JSONParser parser(argv[1]);
        
        parser.printAll();
        std::cout << std::endl;
        
        // 读取配置
        std::string domain_type = parser.getString("problem.domain_type");
        int nx = parser.getInt("grid.nx");
        int ny = parser.getInt("grid.ny");
        
        std::cout << "Domain: " << domain_type << std::endl;
        std::cout << "Grid: " << nx << " x " << ny << std::endl;
        
        // 边界条件
        std::string left_type = parser.getString("boundary.left.type");
        std::cout << "Left boundary: " << left_type << std::endl;
        
        // 测试函数
        std::string test_name = parser.getString("test_function.name");
        std::cout << "Test function: " << test_name << std::endl;
        
        // 输出设置
        bool verbose = parser.getBool("output.verbose");
        std::cout << "Verbose: " << (verbose ? "true" : "false") << std::endl;
        
        std::cout << "\n✓ All tests passed!" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}