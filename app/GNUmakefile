# Makefile by Guillaume Pronost for LEAF example @ 2020/02/09

OSNAME          = $(shell uname -s)

include ../Makefile/Makefile.${OSNAME}

# set compiler options for ROOT
CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += '-fPIC' -std=c++11

INCFLAGS = -I.
INCFLAGS += -I$(WCSIMDIR)/include
INCFLAGS += -I$(BONSAIDIR)/bonsai
INCFLAGS += -I$(LEAFDIR)/leaf

LIBS += ${WCSIMDIR}/libWCSimRoot.so 
LIBS += $(shell root-config --libs --glibs) -lMinuit 
#LIBS_BONSAI += -L${BONSAIDIR} -lWCSimBonsai
#LIBS_BQFIT += -L${LEAFDIR}/lib -lBQFitter
LIBS_BONSAI += ${BONSAIDIR}/libWCSimBonsai.so
LIBS_BQFIT += ${LEAFDIR}/lib/libBQFitter.so


OBJECT = analysis

all: $(OBJECT)

analysis: analysis.o 
	@echo '<< compiling bin analysis >>'
	@$(CXX) -g $(CXXFLAGS) -o $@ $^ $(LIBS) $(LIBS_BONSAI) $(LIBS_BQFIT)

# default rules
.cc.o: 
	@echo '<< compiling' $< '>>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $<

%.o: %.cc %.hh
	@echo '<< compiling' $< '>>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<

%.o: %.cpp
	@echo '<< compiling' $< '>>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<

%.o: %.C %.h
	@echo '<< compiling' $< '>>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<

# cleaner
clean: 
	@echo '<< cleaning >>'
	@rm -f ./*.o
	@rm -f ./lib*.so
	@rm -f ./bin/*
	@rm -f $(OBJECT)
