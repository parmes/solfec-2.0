include Config.mak

# C includes
C_INC=inc/err.h inc/alg.h inc/real.h cpp/version.h

# C++ files
CPP_INC=cpp/solfec.hpp inc/fmt.hpp
CPP_SRC=cpp/tasksys.cpp cpp/solfec.cpp cpp/input.cpp cpp/spline.cpp

# ISPC files
ISPC_SRC=ispc/util.ispc

# ISPC targets
ISPC_TARGETS=sse2,sse4,avx

# MPI C++ compiler
MPICXX=mpicxx

# Python interpeter
PYTHON=python

# Program name
EXE=solfec

# Do the rest
include common.mk
