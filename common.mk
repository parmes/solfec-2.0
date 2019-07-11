ifeq ($(DEBUG),yes)
  CFLAGS=-g -O0 -m64 -fopenmp -DDEBUG
  ISPC=ispc -g -O0 --arch=x86-64 -DDEBUG
else
  CFLAGS=-O2 -m64 -fopenmp
  ISPC=ispc -O2 --arch=x86-64 --woff
endif

ISPC_OBJS4=$(addprefix objs4/, $(ISPC_SRC:.ispc=_ispc.o) $(ISPC_SRC:.ispc=_ispc_sse2.o) $(ISPC_SRC:.ispc=_ispc_sse4.o) $(ISPC_SRC:.ispc=_ispc_avx.o))
ISPC_OBJS8=$(addprefix objs8/, $(ISPC_SRC:.ispc=_ispc.o) $(ISPC_SRC:.ispc=_ispc_sse2.o) $(ISPC_SRC:.ispc=_ispc_sse4.o) $(ISPC_SRC:.ispc=_ispc_avx.o))
ISPC_HEADERS4=$(addprefix objs4/, $(ISPC_SRC:.ispc=_ispc.h))
ISPC_HEADERS8=$(addprefix objs8/, $(ISPC_SRC:.ispc=_ispc.h))
CPP_OBJS4=$(addprefix objs4/, $(CPP_SRC:.cpp=.o))
CPP_OBJS8=$(addprefix objs8/, $(CPP_SRC:.cpp=.o))
C_OBJS4=$(addprefix objs4/, $(C_SRC:.c=.o))
C_OBJS8=$(addprefix objs8/, $(C_SRC:.c=.o))
LIBS=-lm $(PYTHONLIB) $(HDF5LIB)

default: dirs version $(ISPC_HEADERS4) $(EXE)4 $(ISPC_HEADERS8) $(EXE)8

.PHONY: dirs clean print

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

clean:
	/bin/rm -rf objs* *~ $(EXE)4 $(EXE)8 *.dSYM
	find ./ -iname "*.dump" -exec rm '{}' ';'

version:
	python python/version.py

$(EXE)4: $(CPP_OBJS4) $(C_OBJS4) $(ISPC_OBJS4)
	$(MPICXX) $(CFLAGS) -fopenmp -o $@ $^ $(LIBS)

$(EXE)8: $(CPP_OBJS8) $(C_OBJS8) $(ISPC_OBJS8)
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
	$(CXX) -DREALSIZE=4 -Iinc -Iobjs4/ispc $(CFLAGS) $(PYTHONINC) $< -c -o $@

objs8/cpp/input.o: cpp/input.cpp
	$(CXX) -DREALSIZE=8 -Iinc -Iobjs8/ispc $(CFLAGS) $(PYTHONINC) $< -c -o $@

objs4/%.o: %.cpp $(ISPC_HEADERS4)
	$(MPICXX) -DREALSIZE=4 -Iinc -Iobjs4/ispc $(CFLAGS) $< -c -o $@

objs8/%.o: %.cpp $(ISPC_HEADERS8)
	$(MPICXX) -DREALSIZE=8 -Iinc -Iobjs8/ispc $(CFLAGS) $< -c -o $@

objs4/%.o: %.c
	$(MPICC) -DREALSIZE=4 -Iinc -Iobjs4/ispc $(CFLAGS) $< -c -o $@

objs8/%.o: %.c
	$(MPICC) -DREALSIZE=8 -Iinc -Iobjs8/ispc $(CFLAGS) $< -c -o $@
