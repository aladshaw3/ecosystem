#This make file will build c and cpp files with gcc and g++
#Works for any g++ v4.7 or higher

IDIR=include
CC=gcc
CXX=g++
CFLAGS=-I$(IDIR) -Wno-ignored-qualifiers
CXXFLAGS= -I$(IDIR) -Wwrite-strings -Wconversion-null -std=c++11 -fpermissive

ODIR=src/obj
INSDIR=/usr/local/bin

#mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
#current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))

EXE=eco0

_DEPS = config.h lmmin.h lmcurve.h dogfish.h eel.h egret.h error.h \
	finch.h flock.h gsta_opt.h lark.h macaw.h magpie.h mola.h \
	monkfish.h sandbox.h school.h scopsowl_opt.h scopsowl.h shark.h \
	skua_opt.h skua.h yaml_private.h yaml_wrapper.h yaml.h Trajectory.h \
	ui.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = main.o api.o  lmmin.o lmcurve.o dogfish.o dumper.o eel.o egret.o emitter.o \
	error.o finch.o gsta_opt.o lark.o loader.o macaw.o magpie.o \
	mola.o monkfish.o parser.o reader.o sandbox.o scanner.o scopsowl_opt.o \
	scopsowl.o shark.o skua_opt.o skua.o writer.o yaml_wrapper.o Trajectory.o \
	ui.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: src/%.c $(DEPS)
	$(CC) -O3 -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: src/%.cpp $(DEPS)
	$(CXX) -O3 -c -o $@ $< $(CXXFLAGS)

$(EXE): $(OBJ)
	$(CXX) -O3 -o $@ $^ $(CXXFLAGS)

.PHONY: clean install cleanall

clean:
	rm -f $(EXE) $(ODIR)/*.o *~ core $(INCDIR)/*~
install:
	cp $(EXE) $(INSDIR)
	mkdir $(INSDIR)/ecodoc
	cp doc/eco_help_bui.txt $(INSDIR)/ecodoc
	cp doc/eco_help_aui.txt $(INSDIR)/ecodoc
cleanall:
	rm -f $(EXE) $(ODIR)/*.o *~ core $(INCDIR)/*~
	rm -f $(INSDIR)/$(EXE)
	rm -r $(INSDIR)/ecodoc