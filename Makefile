include Config.mak

# C includes
C_INC=inc/err.h inc/alg.h inc/real.h

# C++ files
CPP_INC=cpp/solfec.h cpp/version.h
CPP_SRC=cpp/tasksys.cpp cpp/solfec.cpp cpp/input.cpp cpp/spline.cpp

# ISPC files
ISPC_SRC=ispc/util.ispc

# ISPC targets
ISPC_TARGETS=sse2,sse4,avx

# MPI C++ compiler
MPICXX=mpicxx

# Program name
EXE=solfec

# Do the rest
include common.mk
