#!/usr/bin/env python3
import json
import subprocess
import os
import re

# 备份原始文件
os.system("cp input.json input.json.backup")

# 定义测试组（保持不变）
test_groups = [
    {
        "name": "Group1: Square, All Dirichlet",
        "domain_type": "square",
        "hole": None,
        "boundary": {
            "left": {"type": "dirichlet", "value": "from_test_function"},
            "right": {"type": "dirichlet", "value": "from_test_function"},
            "bottom": {"type": "dirichlet", "value": "from_test_function"},
            "top": {"type": "dirichlet", "value": "from_test_function"}
        }
    },
    {
        "name": "Group2: Square, Mixed Boundary",
        "domain_type": "square",
        "hole": None,
        "boundary": {
            "left": {"type": "dirichlet", "value": "from_test_function"},
            "right": {"type": "dirichlet", "value": "from_test_function"},
            "bottom": {"type": "neumann", "value": "from_test_function"},
            "top": {"type": "neumann", "value": "from_test_function"}
        }
    },
    {
        "name": "Group3: With Hole, All Dirichlet",
        "domain_type": "square_with_hole",
        "hole": {"type": "dirichlet", "value": "from_test_function"},
        "boundary": {
            "left": {"type": "dirichlet", "value": "from_test_function"},
            "right": {"type": "dirichlet", "value": "from_test_function"},
            "bottom": {"type": "dirichlet", "value": "from_test_function"},
            "top": {"type": "dirichlet", "value": "from_test_function"}
        }
    },
    {
        "name": "Group4: With Hole, Square Dirichlet, Hole Neumann",
        "domain_type": "square_with_hole",
        "hole": {"type": "neumann", "value": "from_test_function"},
        "boundary": {
            "left": {"type": "dirichlet", "value": "from_test_function"},
            "right": {"type": "dirichlet", "value": "from_test_function"},
            "bottom": {"type": "dirichlet", "value": "from_test_function"},
            "top": {"type": "dirichlet", "value": "from_test_function"}
        }
    },
    {
        "name": "Group5: With Hole, Square Neumann, Hole Dirichlet",
        "domain_type": "square_with_hole",
        "hole": {"type": "dirichlet", "value": "from_test_function"},
        "boundary": {
            "left": {"type": "neumann", "value": "from_test_function"},
            "right": {"type": "neumann", "value": "from_test_function"},
            "bottom": {"type": "neumann", "value": "from_test_function"},
            "top": {"type": "neumann", "value": "from_test_function"}
        }
    },
    {
        "name": "Group6: With Hole, Square Mixed, Hole Dirichlet",
        "domain_type": "square_with_hole",
        "hole": {"type": "dirichlet", "value": "from_test_function"},
        "boundary": {
            "left": {"type": "dirichlet", "value": "from_test_function"},
            "right": {"type": "dirichlet", "value": "from_test_function"},
            "bottom": {"type": "neumann", "value": "from_test_function"},
            "top": {"type": "neumann", "value": "from_test_function"}
        }
    },
    {
        "name": "Group7: With Hole, Square Mixed, Hole Neumann",
        "domain_type": "square_with_hole",
        "hole": {"type": "neumann", "value": "from_test_function"},
        "boundary": {
            "left": {"type": "dirichlet", "value": "from_test_function"},
            "right": {"type": "dirichlet", "value": "from_test_function"},
            "bottom": {"type": "neumann", "value": "from_test_function"},
            "top": {"type": "neumann", "value": "from_test_function"}
        }
    }
]

# 网格尺寸
grid_sizes = [8, 16,  32, 64]

# 测试函数
test_functions = ["exp_y_sin_x"]

# 圆孔参数
hole_params = {
    "hole_center_x": 0.5,
    "hole_center_y": 0.5,
    "hole_radius": 0.25
}

# 读取基础配置
with open("input.json", "r") as f:
    base_config = json.load(f)

# 打开结果文件
with open("test_results.txt", "w") as results_file:
    results_file.write("=== Poisson Solver Test Results ===\n\n")
    
    for test_func in test_functions:
        for group in test_groups:
            for nx in grid_sizes:
                print(f"\nTesting: {group['name']}, {test_func}, n={nx}")
                
                results_file.write(f"\n{group['name']}\n")
                results_file.write(f"  Test function: {test_func}, Grid: {nx} x {nx}\n")
                
                # 修改配置
                base_config["problem"]["domain_type"] = group["domain_type"]
                
                if group["domain_type"] == "square_with_hole":
                    base_config["problem"]["hole_center_x"] = hole_params["hole_center_x"]
                    base_config["problem"]["hole_center_y"] = hole_params["hole_center_y"]
                    base_config["problem"]["hole_radius"] = hole_params["hole_radius"]
                    base_config["boundary"]["hole"] = group["hole"]
                else:
                    if "hole" in base_config["boundary"]:
                        del base_config["boundary"]["hole"]
                
                base_config["boundary"]["left"] = group["boundary"]["left"]
                base_config["boundary"]["right"] = group["boundary"]["right"]
                base_config["boundary"]["bottom"] = group["boundary"]["bottom"]
                base_config["boundary"]["top"] = group["boundary"]["top"]
                
                base_config["grid"]["n"] = nx
                
                base_config["test_function"]["name"] = test_func
                base_config["output"]["verbose"] = False
                base_config["output"]["save_solution"] = False
                
                # 保存配置
                with open("input_temp.json", "w") as f:
                    json.dump(base_config, f, indent=4)
                
                # 运行求解器
                result = subprocess.run(["./solver", "input_temp.json"], 
                                       capture_output=True, text=True)
                
                # 提取误差范数
                output = result.stdout
                
                l1_match = re.search(r'L1 norm:\s+([\d\.eE+-]+)', output)
                l2_match = re.search(r'L2 norm:\s+([\d\.eE+-]+)', output)
                linf_match = re.search(r'L∞ norm:\s+([\d\.eE+-]+)', output)
                
                if l1_match and l2_match and linf_match:
                    l1 = l1_match.group(1)
                    l2 = l2_match.group(1)
                    linf = linf_match.group(1)
                    result_line = f"  L1 norm: {l1}, L2 norm: {l2}, L∞ norm: {linf}\n"
                    print(result_line.strip())
                    results_file.write(result_line)
                else:
                    print("  Error: Could not extract error norms")
                    results_file.write("  Error: Could not extract error norms\n")
                
                # 删除临时文件
                os.remove("input_temp.json")

# 恢复原始文件
os.system("mv input.json.backup input.json")

print("\n" + "="*60)
print("All tests completed! Results saved to test_results.txt")
print("="*60)
