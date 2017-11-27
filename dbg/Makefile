#This make file will build c and cpp files with gcc and g++
#Works for any g++ v4.7 or higher

IDIR=../include
CC=gcc
CXX=g++
CFLAGS=-I$(IDIR) -Wno-discarded-qualifiers
CXXFLAGS= -I$(IDIR) -Wwrite-strings -Wconversion-null -std=c++11 -fpermissive

ODIR=obj

DBG=eco-dbg

_DEPS = config.h lmmin.h lmcurve.h dogfish.h eel.h egret.h error.h \
	finch.h flock.h gsta_opt.h lark.h macaw.h magpie.h mola.h \
	monkfish.h sandbox.h school.h scopsowl_opt.h scopsowl.h shark.h \
	skua_opt.h skua.h yaml_private.h yaml_wrapper.h yaml.h Trajectory.h \
	ui.h dove.h crow.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o api.o  lmmin.o lmcurve.o dogfish.o dumper.o eel.o egret.o emitter.o \
	error.o finch.o gsta_opt.o lark.o loader.o macaw.o magpie.o \
	mola.o monkfish.o parser.o reader.o sandbox.o scanner.o scopsowl_opt.o \
	scopsowl.o shark.o skua_opt.o skua.o writer.o yaml_wrapper.o Trajectory.o \
	ui.o dove.o crow.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: ../src/%.c $(DEPS)
	$(CC) -g3 -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: ../src/%.cpp $(DEPS)
	$(CXX) -g3 -c -o $@ $< $(CXXFLAGS)

$(DBG): $(OBJ)
	$(CXX) -g3 -o $@ $^ $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(DBG) $(ODIR)/*.o *~ core $(INCDIR)/*~