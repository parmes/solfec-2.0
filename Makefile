include Config.mak

# C++ files
CPP_SRC=cpp/tasksys.cpp cpp/solfec.cpp cpp/input.cpp

# ISPC files
ISPC_SRC=ispc/simd_util.ispc

# ISPC targets
ISPC_TARGETS=sse2,sse4,avx

# MPI C++ compiler
MPICXX=mpicxx

# Program name
EXE=solfec

# Do the rest
include common.mk
