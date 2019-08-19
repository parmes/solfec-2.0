include Config.mak

# C includes
C_INC=inc/err.h inc/alg.h inc/real.h cpp/version.h

# C++ files
CPP_INC=cpp/solfec.hpp inc/fmt.hpp cpp/compute.hpp cpp/mesh.hpp cpp/part.hpp\
        cpp/mutex.hpp cpp/lock.hpp cpp/ga.hpp
CPP_SRC=cpp/tasksys.cpp cpp/solfec.cpp cpp/input.cpp cpp/spline.cpp cpp/output.cpp\
        cpp/compute.cpp cpp/mesh.cpp cpp/part.cpp cpp/mutex.cpp cpp/lock.cpp cpp/ga.cpp\
	cpp/dynlb.cpp

# ISPC files
ISPC_SRC=ispc/util.ispc ispc/alloc.ispc ispc/part.ispc ispc/rcb.ispc

# ISPC targets
ISPC_TARGETS=sse2,sse4,avx

# MPI C++ compiler
MPICXX=mpicxx -std=c++17

# Python interpeter
PYTHON=python

# Program name
EXE=solfec

# Do the rest
include common.mk
