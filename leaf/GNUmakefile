
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

SOURCES 	= LEAF.cc
HEADERS 	= LEAF.hh

include ../Makefile/Makefile.${OSNAME}

ROOT_VERSION = $(shell root-config --version | tr -d '/' )


# set compiler options for ROOT
CXXFLAGS_BASE += $(shell root-config --cflags)
CXXFLAGS_BASE += '-fPIC' -Wall -Wpedantic -Wno-long-long -DDATAMODEL_LITE

CXXFLAGS = $(CXXFLAGS_BASE)
#CXXFLAGS_ROOT = $(CXXFLAGS_BASE)

ifeq ($(shell echo "$(ROOT_VERSION) >= 6.0" | bc -l ), 1)
MAKE_LEAF_WITH_ROOT6 = 1
CXXFLAGS += -std=c++11 -DROOT6
#CXXFLAGS_ROOT += -std=c++11 -DROOT6
else
CXXFLAGS += -std=c++11 -DROOT5
#CXXFLAGS_ROOT += -DROOT5
endif

INCFLAGS = -I. -I$(shell root-config --incdir)
INCFLAGS += -I$(WCSIMDIR)/include
INCFLAGS += -I./DataModel

LIBS += -L${WCSIMDIR} -lWCSimRoot
LIBS += $(shell root-config --libs) -lMinuit 

LIB_DATAMODEL = ../lib/libDataModelLite.so
LIB_MANAGER = ../lib/libHKManager.so
LIB_FITTER = ../lib/libLEAF.so
LIB_ASTROANALYSIS = ../lib/libHKAstroAnalysis.so

LIBS_OBJECTS = $(LIB_DATAMODEL) $(LIB_MANAGER) $(LIB_FITTER) #$(LIB_ASTROANALYSIS)


DATAMODEL_SOURCES = $(wildcard DataModel/*.cpp)
DATAMODEL_SRC_NOP = $(DATAMODEL_SOURCES:DataModel/%=./%)

DATAMODEL_HEADERS = $(DATAMODEL_SOURCES:%.cpp=%.h)
DATAMODEL_OBJECTS = $(DATAMODEL_SRC_NOP:%.cpp=%.o) #DataModelRootDict.o


all: .depend ROOT_CHECK $(LIBS_OBJECTS)

ROOT_CHECK:
ifeq ($(MAKE_LEAF_WITH_ROOT6),1)
	@echo '<< compile with ROOT6 >>'
else
	@echo '<< compile with ROOT5 >>'
endif

# library
$(LIB_FITTER): LEAF.o 
	@if [ ! -d ../lib ]; then \
		mkdir ../lib; \
	fi
	@echo '<< generating ' $@ ' >>'
	@$(SOMAKER) $(SOFLAGS) -o $@ $^ $(LIBS)
	
$(LIB_MANAGER): HKManager.o 
	@echo '<< generating ' $@ ' >>'
	@if [ ! -d ../lib ]; then \
		mkdir ../lib; \
	fi
	@$(SOMAKER) $(SOFLAGS) -o $@ $^ $(LIBS)
	
$(LIB_ASTROANALYSIS): HKAstroAnalysis.o 
	@echo '<< generating ' $@ ' >>'
	@if [ ! -d ../lib ]; then \
		mkdir ../lib; \
	fi
	@$(SOMAKER) $(SOFLAGS) -o $@ $^ $(LIBS)
	
#$(LIB_DATAMODEL): DataModelRootDict.cpp $(DATAMODEL_OBJECTS) 
$(LIB_DATAMODEL): $(DATAMODEL_OBJECTS) 
	@echo '<< generating ' $@ ' >>'
	@$(SOMAKER) $(SOFLAGS) -o $@  $(DATAMODEL_OBJECTS)  $(LIBS)
	
# Not needed anymore
#DataModelRootDict.cpp:
#	@echo '<< generating ' $@ ' >>'
#	@rootcint -f DataModelRootDict.cpp -c $(CXXFLAGS_ROOT) -p $(DATAMODEL_HEADERS) $(PWD)/DataModel/LinkDef.h
#ifeq ($(MAKE_LEAF_WITH_ROOT6),1)
#	@echo '<< install DataModelRootDict_rdict.pcm >>'
#	@cp DataModelRootDict_rdict.pcm ../lib/.
#endif
	
# DataModel objects
%.o: DataModel/%.cpp DataModel/%.h
	@echo '<< compiling DataModel object: ' $@ ' >>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<

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


###################################################
# auto-generate dependencies using g++ -M
.depend: $(SOURCES) $(HEADERS)
	@#rm .depend
	@echo '<< checking dependencies >>'
	@g++ -M $(CXXFLAGS) $(INCFLAGS) $(SOURCES) > .depend

-include .depend

# cleaner
clean: 
	@echo '<< cleaning >>'
	@rm -f .depend
	@rm -f ./*.o
	@rm -f ./../lib/lib*.so
	@rm -f ./../lib/*Dict*
	@rm -f ./*Dict*
	
