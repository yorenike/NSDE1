# 编译器设置
CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2

# 源文件
SOURCES = JSONParser.cpp Grid.cpp TestFunction.cpp FDDiscretization.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# 可执行文件
TARGET = solver

# 默认目标
all: $(TARGET)

# 主求解器
$(TARGET): main.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ main.o $(OBJECTS)

# 编译规则
main.o: main.cpp JSONParser.h Grid.h TestFunction.h FDDiscretization.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

JSONParser.o: JSONParser.cpp JSONParser.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

Grid.o: Grid.cpp Grid.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

TestFunction.o: TestFunction.cpp TestFunction.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

FDDiscretization.o: FDDiscretization.cpp FDDiscretization.h Grid.h TestFunction.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 运行主程序
run: $(TARGET)
	./$(TARGET) input.json

# 清理
clean:
	rm -f *.o $(TARGET)

# 帮助
help:
	@echo "Available targets:"
	@echo "  make          - 编译程序"
	@echo "  make run      - 编译并运行主程序"
	@echo "  make clean    - 清理目标文件和可执行文件"
	@echo "  make help     - 显示帮助信息"

.PHONY: all run clean help