MAIN = fitParam

# simple Makefile to build the JunoARCon.exe
#ROOTLIBS   = $(shell ~/tools/bin/root-config --libs)
#ROOTGLIBS   = $(shell ~/tools/bin/root-config --glibs)
#ROOTCFLAGS = $(shell ~/tools/bin/root-config --cflags)
ROOTLIBS   = $(shell  root-config --libs)
ROOTGLIBS   = $(shell root-config --glibs)
ROOTCFLAGS = $(shell  root-config --cflags)

#work with linux
#CXX       = clang++ 

# options for compiler
CXXFLAGS  = -Wall
CXXFLAGS  += -I. -I./include
CXXFLAGS  += $(ROOTCFLAGS)

# options for linker
LDFLAGS  = -O2 -Wall -fPIC
LIBS      = $(ROOTLIBS) 
LIBS += -lMinuit

# generate JunoARCon.exe
SRCS = $(wildcard src/*.cc)

HHS = $(SRCS:src/%.cc=include/%.hh) #$(wildcard include/*.hh)

#SRCS += JunoARCon.cc

OBJS = $(SRCS:src/%.cc=obj/%.o)

obj/%.o : src/%.cc
	$(CXX) -o $@ $(CXXFLAGS) -c $<

all : $(MAIN).exe

$(MAIN).o : $(MAIN).cc $(HHS)
	$(CXX) -o $@ $(CXXFLAGS) -I./include -c $<

$(OBJS) : obj/%.o:src/%.cc include/%.hh

$(MAIN).exe : %.exe:%.o ${OBJS} 
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

clean:
		@rm -f ${OBJS} 
		@rm -f *.o 
help:
	@echo "SRCS=$(SRCS)"
	@echo "OBJS=$(OBJS)"
	@echo "HHS=$(HHS)"
