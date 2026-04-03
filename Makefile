# 编译器设置
CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2

# 目录设置
SRC_DIR = src
BUILD_DIR = build
RESULTS_DIR = results

# 源文件（在 src 目录下）
SOURCES = $(SRC_DIR)/JSONParser.cpp \
          $(SRC_DIR)/Grid.cpp \
          $(SRC_DIR)/TestFunction.cpp \
          $(SRC_DIR)/FDDiscretization.cpp \
          $(SRC_DIR)/LinearSolver.cpp

# 目标文件（放在 build 目录下）
OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# 可执行文件（放在根目录）
TARGET = solver

# 默认目标
all: $(BUILD_DIR) $(TARGET)

# 创建 build 目录
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# 创建 results 目录
$(RESULTS_DIR):
	mkdir -p $(RESULTS_DIR)

# 主求解器
$(TARGET): $(BUILD_DIR)/main.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# 编译 main.o（依赖 build 目录）
$(BUILD_DIR)/main.o: $(SRC_DIR)/main.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 编译各个模块（依赖 build 目录）
$(BUILD_DIR)/JSONParser.o: $(SRC_DIR)/JSONParser.cpp $(SRC_DIR)/JSONParser.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/Grid.o: $(SRC_DIR)/Grid.cpp $(SRC_DIR)/Grid.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/TestFunction.o: $(SRC_DIR)/TestFunction.cpp $(SRC_DIR)/TestFunction.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/FDDiscretization.o: $(SRC_DIR)/FDDiscretization.cpp $(SRC_DIR)/FDDiscretization.h $(SRC_DIR)/Grid.h $(SRC_DIR)/TestFunction.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/LinearSolver.o: $(SRC_DIR)/LinearSolver.cpp $(SRC_DIR)/LinearSolver.h | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 运行主程序（自动创建 results 目录）
run: $(TARGET) $(RESULTS_DIR)
	./$(TARGET) input.json

# 清理
clean:
	rm -rf $(BUILD_DIR) $(TARGET)

# 完全清理（包括结果）
cleanall: clean
	rm -rf $(RESULTS_DIR)

# 帮助
help:
	@echo "Available targets:"
	@echo "  make          - 编译程序"
	@echo "  make run      - 编译并运行主程序"
	@echo "  make clean    - 清理目标文件和可执行文件"
	@echo "  make cleanall - 清理所有（包括结果文件）"
	@echo "  make help     - 显示帮助信息"

.PHONY: all run clean cleanall help