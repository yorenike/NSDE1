CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2
LAPACK_LIBS = -llapack -lblas

# 源文件
SOURCES = JSONParser.cpp Grid.cpp TestFunction.cpp FDDiscretization.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# 可执行文件
TARGET = solver
CONV_TEST = convergence_test   # 改成 convergence_test

# 默认目标
all: $(TARGET) $(CONV_TEST)

# 主求解器
$(TARGET): main.o $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ main.o $(OBJECTS) $(LAPACK_LIBS)

# 收敛性测试
$(CONV_TEST): convergence_test.o $(OBJECTS)   # 改成 convergence_test.o
	$(CXX) $(CXXFLAGS) -o $@ convergence_test.o $(OBJECTS) $(LAPACK_LIBS)

# 编译规则
%.o: %.cpp JSONParser.h Grid.h TestFunction.h FDDiscretization.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 运行主程序
run: $(TARGET)
	./$(TARGET) input.json

# 运行收敛性测试
test: $(CONV_TEST)
	./$(CONV_TEST)

# 清理
clean:
	rm -f *.o $(TARGET) $(CONV_TEST)

# 帮助
help:
	@echo "Available targets:"
	@echo "  make          - 编译主程序和测试程序"
	@echo "  make run      - 编译并运行主程序"
	@echo "  make test     - 编译并运行收敛性测试"
	@echo "  make clean    - 清理目标文件和可执行文件"