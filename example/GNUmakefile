# Makefile by Guillaume Pronost for LEAF example @ 2020/02/09

OSNAME          = $(shell uname -s)

include ../Makefile/Makefile.${OSNAME}


# set compiler options for ROOT
CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += -fPIC -std=c++17 -Wall -Wpedantic -Wno-long-long

INCFLAGS = -I. 
#INCFLAGS += -I$(shell root-config --incdir) #included in --cflags
INCFLAGS += -I$(LEAFDIR)/libWCSIM/include
INCFLAGS += -I$(LEAFDIR)/leaf
INCFLAGS += -I$(LEAFDIR)/leaf/DataModel-lite

# LIBS += -L${WCSIMDIR} -lWCSimRoot 
LIBS += -L$(LEAFDIR)/libWCSIM -lWCSimRoot
LIBS += $(shell root-config --libs) -lMinuit 
# LIBS_BONSAI += -L${BONSAIDIR} -lWCSimBonsai
LIBS_LEAF += -L${LEAFDIR}/lib -lDataModelLite -lHKManager -lLEAF #-lHKAstroAnalysis 


OBJECT = analysis

all: $(OBJECT)
	
analysis: analysis.o 
	@echo '<< compiling bin analysis >>'
	@$(CXX) -g $(CXXFLAGS) -o $@ $^ $(LIBS) $(LIBS_LEAF)
		
# default rules
.cc.o: 
	@echo '<< compiling' $< '>>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $< 

%.o: %.cpp %.hpp
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
	
