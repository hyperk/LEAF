# Makefile by Guillaume Pronost for LEAF example @ 2020/02/09

OSNAME          = $(shell uname -s)
LEAFDIR         = ..

include ../Makefile/Makefile.${OSNAME}

# set compiler options for ROOT
CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += -fPIC -std=c++17 -Wall -Wpedantic -Wno-long-long
ROOT_VERSION = $(shell root-config --version | cut -d'.' -f1,2)
ifeq ($(shell echo "$(ROOT_VERSION) >= 6.0" | bc -l), 1)
	HK_USE_ROOT7 = 1
	CXXFLAGS += -DHK_USE_ROOT7
endif

INCFLAGS = -I. 
#INCFLAGS += -I$(shell root-config --incdir) #included in --cflags
INCFLAGS += -I$(LEAFDIR)/libWCSIM/include
INCFLAGS += -I$(LEAFDIR)/leaf
INCFLAGS += -I$(LEAFDIR)/leaf/DataModel-lite

LIBS += -L$(LEAFDIR)/libWCSIM -lWCSimRoot
LIBS += $(shell root-config --libs) -lMinuit 
# LIBS_BONSAI += -L${BONSAIDIR} -lWCSimBonsai
LIBS_LEAF += -L${LEAFDIR}/lib -lLEAF -lHKManager -lDataModelLite

OBJECT = analysis

all: ROOT_CHECK $(OBJECT)

ROOT_CHECK:
ifeq ($(HK_USE_ROOT7),1)
	@echo '<< compile with ROOT7 >>'
else
	@echo '<< compile with old ROOT >>'
	@echo $(DATAMODEL_SOURCES)	
endif
	
analysis: analysis.o 
	@echo '<< compiling bin $@ >>'
	@$(CXX) -g $(CXXFLAGS) $^ $(LIBS) $(LIBS_LEAF) -o $@
		
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
	
