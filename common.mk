ifeq ($(DEBUG),yes)
  CFLAGS=-g -O0 -m64 -fopenmp -DDEBUG
  ISPC=ispc -g -O0 --arch=x86-64 -DDEBUG
  OUTS = $(ISPC_HEADERS4) $(EXE)4
else
  CFLAGS=-O2 -m64 -fopenmp
  ISPC=ispc -O2 --arch=x86-64 --woff
  OUTS = $(ISPC_HEADERS4) $(EXE)4 $(ISPC_HEADERS8) $(EXE)8
endif

ISPC_OBJS4=$(addprefix objs4/, $(ISPC_SRC:.ispc=_ispc.o) $(ISPC_SRC:.ispc=_ispc_sse2.o) $(ISPC_SRC:.ispc=_ispc_sse4.o) $(ISPC_SRC:.ispc=_ispc_avx.o))
ISPC_OBJS8=$(addprefix objs8/, $(ISPC_SRC:.ispc=_ispc.o) $(ISPC_SRC:.ispc=_ispc_sse2.o) $(ISPC_SRC:.ispc=_ispc_sse4.o) $(ISPC_SRC:.ispc=_ispc_avx.o))
ISPC_HEADERS4=$(addprefix objs4/, $(ISPC_SRC:.ispc=_ispc.h))
ISPC_HEADERS8=$(addprefix objs8/, $(ISPC_SRC:.ispc=_ispc.h))
CPP_OBJS4=$(addprefix objs4/, $(CPP_SRC:.cpp=.o))
CPP_OBJS8=$(addprefix objs8/, $(CPP_SRC:.cpp=.o))
C_OBJS4=$(addprefix objs4/, $(C_SRC:.c=.o))
C_OBJS8=$(addprefix objs8/, $(C_SRC:.c=.o))
LIBS=-lm $(PYTHONLIB) $(HDF5LIB) -Llib/metis -lgklib -lmetis

default: libmetis inc/taskflow inc/amgcl dirs version $(OUTS)

.PHONY: dirs clean print test

libmetis: lib/metis
	cd lib/metis && make

lib/metis:
	@echo "Unzipping local copy of METIS graph partitioning library..."
	cd lib && unzip metis.zip
	@echo "METIS extracted."

inc/taskflow:
	@echo "Downloading taskflow headers from GitHub into inc/taskflow..."
	cd inc && make taskflow
	@echo "You can update inc/taskflow by issuing 'make taskflow' within inc/"

inc/amgcl:
	@echo "Downloading amgcl headers from GitHub into inc/amgcl..."
	cd inc && make amgcl
	@echo "You can update inc/amgcl by issuing 'make amgcl' within inc/"

test:
	cd tests && $(PYTHON) -m unittest discover

print:
	@echo $(ISPC_HEADERS4)
	@echo $(CPP_OBJS4)
	@echo $(C_OBJS4)
	@echo $(ISPC_OBJS4)
	@echo $(ISPC_HEADERS8)
	@echo $(CPP_OBJS8)
	@echo $(C_OBJS8)
	@echo $(ISPC_OBJS8)

dirs:
	/bin/mkdir -p objs4/
	/bin/mkdir -p objs8/
	/bin/mkdir -p objs4/cpp
	/bin/mkdir -p objs8/cpp
	/bin/mkdir -p objs4/ispc
	/bin/mkdir -p objs8/ispc

del:
	find ./ -iname "*.dump" -exec rm '{}' ';'
	find ./ -iname "debug_*.txt" -exec rm '{}' ';'
	rm -fr tests/*.h5
	rm -fr tests/*.xmf

clean: del
	/bin/rm -rf objs* *~ $(EXE)4 $(EXE)8 *.dSYM

cleanall: clean
	cd lib/metis && make clean

version:
	python python/version.py

$(EXE)4: $(ISPC_OBJS4) $(CPP_OBJS4) $(C_OBJS4)
	$(MPICXX) $(CFLAGS) -fopenmp -o $@ $^ $(LIBS)

$(EXE)8: $(ISPC_OBJS8) $(CPP_OBJS8) $(C_OBJS8)
	$(MPICXX) $(CFLAGS) -fopenmp -o $@ $^ $(LIBS)

objs4/%_ispc.h objs4/%_ispc.o objs4/%_ispc_sse2.o objs4/%_ispc_sse4.o objs4/%_ispc_avx.o: %.ispc
	$(ISPC) -DREALSIZE=4 --target=$(ISPC_TARGETS) -Iinc $< -o objs4/$*_ispc.o -h objs4/$*_ispc.h

objs8/%_ispc.h objs8/%_ispc.o objs8/%_ispc_sse2.o objs8/%_ispc_sse4.o objs8/%_ispc_avx.o: %.ispc
	$(ISPC) -DREALSIZE=8 --target=$(ISPC_TARGETS) -Iinc $< -o objs8/$*_ispc.o -h objs8/$*_ispc.h

objs4/cpp/tasksys.o: cpp/tasksys.cpp
	$(MPICXX) $(CFLAGS) -D ISPC_USE_OMP $< -c -o $@

objs8/cpp/tasksys.o: cpp/tasksys.cpp
	$(MPICXX) $(CFLAGS) -D ISPC_USE_OMP $< -c -o $@

objs4/cpp/input.o: cpp/input.cpp
	$(MPICXX) -DREALSIZE=4 -Iinc -Iobjs4/ispc $(CFLAGS) $(PYTHONINC) $< -c -o $@

objs8/cpp/input.o: cpp/input.cpp
	$(MPICXX) -DREALSIZE=8 -Iinc -Iobjs8/ispc $(CFLAGS) $(PYTHONINC) $< -c -o $@

objs4/cpp/output.o: cpp/output.cpp
	$(MPICXX) -DREALSIZE=4 -Iinc -Iobjs4/ispc $(CFLAGS) $(HDF5INC) $< -c -o $@

objs8/cpp/output.o: cpp/output.cpp
	$(MPICXX) -DREALSIZE=8 -Iinc -Iobjs8/ispc $(CFLAGS) $(HDF5INC) $< -c -o $@

objs4/%.o: %.cpp $(ISPC_HEADERS4) $(CPP_INC) $(C_INC)
	$(MPICXX) -DREALSIZE=4 -Iinc -Ilib -Iobjs4/ispc $(CFLAGS) $< -c -o $@

objs8/%.o: %.cpp $(ISPC_HEADERS8) $(CPP_INC) $(C_INC)
	$(MPICXX) -DREALSIZE=8 -Iinc -Ilib -Iobjs8/ispc $(CFLAGS) $< -c -o $@

objs4/%.o: %.c $(C_INC)
	$(MPICC) -DREALSIZE=4 -Iinc -Iobjs4/ispc $(CFLAGS) $< -c -o $@

objs8/%.o: %.c $(C_INC)
	$(MPICC) -DREALSIZE=8 -Iinc -Iobjs8/ispc $(CFLAGS) $< -c -o $@
