DIR = $(shell pwd)

DEL           = Delphes
#CXX           = ~/.vim/bundle/ClangComplete/bin/cc_args.py g++
CXX           = g++
#CXXFLAGS      = -Wno-deprecated $(shell root-config --cflags) -I../$(DEL)/external/ -I../$(DEL) -std=c++0x 
CCFLAGS       = $(shell root-config --cflags) -I$(DEL)/external/ -I$(DEL) -std=c++0x
CXXFLAGS      = $(shell root-config --cflags) -I$(DEL)/external/ -I$(DEL)
LDFLAGS       = 
LD            = g++

LIBS = $(shell root-config --glibs)
LIBS += -LDelphes -lDelphes

OBJS          = DelWeight.o main.o
PROGRAM       = main

%.o : %.C $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

%.o : %.cc $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CCFLAGS) -c $<

%.o : %.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $<

$(PROGRAM): $(OBJS)
	@echo "Linking $(PROGRAM) ..."
	@$(LD) $(OBJS) $(LIBS) -o $(PROGRAM)
	@echo "done"
