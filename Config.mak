# ISPC (http://ispc.github.io) and MPICXX compilers are assumed to be in the PATH

# Python paths (used to parse input files)
PYTHONINC=-I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
PYTHONLIB=-L/opt/local/lib -lpython2.7

# HDF5 paths (used to write output files)
HDF5INC=-I/opt/local/include
HDF5LIB=-L/opt/local/lib -lhdf5 -lhdf5_hl

# Debug version
DEBUG=yes
