INTELROOT = /opt/intel
MKLROOT   = $(INTELROOT)/mkl
EIGENDIR  = /usr/local/include/eigen3

CXX      = g++-11
MEX      = mex
MWRAP    = mwrap
LIBS     = -lgomp -lmkl_intel_lp64 -lmkl_core -lm -ldl -lmkl_sequential
CXXFLAGS = -fopenmp -std=c++20 -O3 -DNDEBUG -march=native -fPIC -fno-math-errno
CXXFLAGS += -DMKL_THREADING_LAYER=sequential -DMKL_NUM_THREADS=1 -DMKL_DYNAMIC="FALSE" -DMKL_DOMAIN_NUM_THREADS="MKL_DOMAIN_ALL=1, MKL_DOMAIN_BLAS=1" -DOMP_NUM_THREADS=1 -DOMP_DYNAMIC="FALSE"
LDFLAGS  = -bundle
MEXFLAGS = -I$(EIGENDIR) -I$(MKLROOT)/include -largeArrayDims -DR2008OO CXX=$(CXX) CXXFLAGS='$(CXXFLAGS)' CXXOPTIMFLAGS='' LDFLAGS='$(LDFLAGS)'

# MKL directories on Mac and Linux are different
OS := $(shell uname)
ifeq ($(OS),Darwin)
    LIBDIR = -L$(INTELROOT)/lib -L$(MKLROOT)/lib
    LDFLAGS += -Wl,-rpath,$(MKLROOT)/lib,-rpath,$(INTELROOT)/lib
else ifeq ($(OS),Linux)
    LIBDIR = -L$(INTELROOT)/lib/intel64 -L$(MKLROOT)/lib/intel64 -Wl,--no-as-needed
endif
