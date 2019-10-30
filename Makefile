include Config.mak

# C includes
C_INC=inc/err.h inc/alg.h inc/real.h cpp/version.h

# C++ files
CPP_INC=cpp/solfec.hpp inc/fmt.hpp cpp/compute.hpp cpp/mesh.hpp cpp/part.hpp\
        cpp/mutex.hpp cpp/lock.hpp cpp/ga.hpp cpp/debug.hpp cpp/timestep.hpp
CPP_SRC=cpp/tasksys.cpp cpp/solfec.cpp cpp/input.cpp cpp/spline.cpp cpp/output.cpp\
        cpp/compute.cpp cpp/mesh.cpp cpp/part.cpp cpp/mutex.cpp cpp/lock.cpp cpp/ga.cpp\
	cpp/dynlb.cpp cpp/debug.cpp cpp/timestep.cpp

# ISPC files
ISPC_SRC=ispc/util.ispc ispc/alloc.ispc ispc/part.ispc ispc/rcb.ispc\
         ispc/hex0.ispc ispc/wed0.ispc ispc/pyr0.ispc ispc/tet0.ispc

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
