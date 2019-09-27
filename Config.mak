# ISPC (http://ispc.github.io) and MPICXX compilers are assumed to be in the PATH

# Python {2 or 3} paths (Python is used to parse input files)
#PYTHONINC=-I/opt/local/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
#PYTHONLIB=-L/opt/local/lib -lpython2.7
PYTHONINC=-I/opt/local/Library/Frameworks/Python.framework/Versions/3.7/include/python3.7m
PYTHONLIB=-L/opt/local/Library/Frameworks/Python.framework/Versions/3.7/lib/python3.7/config-3.7m-darwin -lpython3.7m -ldl -framework CoreFoundation

# HDF5 paths (HDF5 is used to write output files)
HDF5INC=-I/usr/local/include
HDF5LIB=-L/usr/local/lib -lhdf5 -lhdf5_hl

# Debug version
DEBUG=yes
