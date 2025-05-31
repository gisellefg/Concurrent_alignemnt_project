# Compiler and flags
CXX = g++
NVCC = nvcc

CXXFLAGS = -std=c++17 -O3
NVCCFLAGS = -std=c++17 -O3

# Target executable
TARGET = align_benchmark

# Source files
CPU_SRC = gotohCPU.cpp
CUDA_SRC = gotohCUDA.cu
MAIN_SRC = main.cpp

# Object files
CPU_OBJ = gotohCPU.o
CUDA_OBJ = gotohCUDA.o
MAIN_OBJ = main.o

all: $(TARGET)

# Compile CPU source
$(CPU_OBJ): $(CPU_SRC) gotohCPU.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile CUDA source
$(CUDA_OBJ): $(CUDA_SRC) gotohCUDA.h
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# Compile main source
$(MAIN_OBJ): $(MAIN_SRC) gotohCPU.h gotohCUDA.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link all
$(TARGET): $(MAIN_OBJ) $(CPU_OBJ) $(CUDA_OBJ)
	$(NVCC) $(NVCCFLAGS) $^ -o $@

clean:
	rm -f $(TARGET) *.o

.PHONY: all clean
