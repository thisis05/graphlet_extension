LLVM_DIR = /opt/homebrew/opt/llvm
OPENMP_LIB ?= 1

all: base_4_d_motif random_4_d_motif 4_d_motif dTri.exe

CXXFLAGS = -O3
ifeq ($(OPENMP_LIB),1)
    CXXFLAGS += -fopenmp
endif
CXXFLAGS += -ffast-math -funroll-loops -fno-strict-aliasing \
          -fomit-frame-pointer -march=native

CXX = $(LLVM_DIR)/bin/clang++

base_4_d_motif: base_4_d_motif.cpp
	$(CXX) $(CXXFLAGS) -o base_4_d_motif base_4_d_motif.cpp

random_4_d_motif: random_4_d_motif.cpp
	$(CXX) $(CXXFLAGS) -o random_4_d_motif random_4_d_motif.cpp

4_d_motif: 4_d_motif.cpp
	$(CXX) $(CXXFLAGS) -o 4_d_motif 4_d_motif.cpp

dTri.exe: dTri.cpp vertex.cpp edge.cpp vertex.h edge.h 
	$(CXX) $(CXXFLAGS) -o dTri.exe dTri.cpp vertex.cpp edge.cpp

clean:
	@echo "Removing executables"
	rm -rf base_4_d_motif random_4_d_motif 4_d_motif dTri.exe
