
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

SOURCES 	= LEAF.cpp
HEADERS 	= LEAF.hpp

include ../Makefile/Makefile.${OSNAME}

# set compiler options for ROOT
CXXFLAGS_BASE += $(shell root-config --cflags)
CXXFLAGS_BASE += '-fPIC' -Wall -Wpedantic -Wno-long-long -DDATAMODEL_LITE

CXXFLAGS = $(CXXFLAGS_BASE) -std=c++17

ROOT_VERSION = $(shell root-config --version | cut -d'.' -f1,2)
ifeq ($(shell echo "$(ROOT_VERSION) >= 6.0" | bc -l), 1)
	HK_USE_ROOT7 = 1
	CXXFLAGS += -DHK_USE_ROOT7
endif

INCFLAGS = -I. -I$(shell root-config --incdir)
INCFLAGS += -I../libWCSIM/include
INCFLAGS += -I./DataModel-lite

LIBS += -L../libWCSIM -lWCSimRoot
LIBS += $(shell root-config --libs) -lMinuit 

LIB_DATAMODEL = ../lib/libDataModelLite.so
LIB_MANAGER = ../lib/libHKManager.so
LIB_FITTER = ../lib/libLEAF.so
PCM_LIB_DATAMODEL = ../lib/libDataModelLite_rdict.pcm
PCM_DATAMODEL = ../lib/DataModelRootDict_rdict.pcm

LIBS_OBJECTS = $(LIB_DATAMODEL) $(LIB_MANAGER) $(LIB_FITTER) 

DATAMODEL_SOURCES = $(filter-out DataModel-lite/DataModelRootDict.cpp, $(wildcard DataModel-lite/*.cpp))
DATAMODEL_SRC_NOP = $(DATAMODEL_SOURCES:DataModel-lite/%=./%)

DATAMODEL_HEADERS = $(DATAMODEL_SOURCES:%.cpp=%.hpp)
DATAMODEL_OBJECTS = $(DATAMODEL_SRC_NOP:%.cpp=%.o) DataModelRootDict.o

all: .depend ROOT_CHECK DIRECTORY_CHECK $(LIBS_OBJECTS) $(PCM_DATAMODEL) 

ROOT_CHECK:
ifeq ($(HK_USE_ROOT7),1)
	@echo '<< compile with ROOT7 >>'
else
	@echo '<< compile with old ROOT >>'
	@echo $(DATAMODEL_SOURCES)	
endif

DIRECTORY_CHECK:
	@if [ ! -d ../lib ]; then \
		mkdir ../lib; \
	fi

# library
$(LIB_FITTER): LEAF.o LEAFUtilities.o SplineUtilities.o
	@echo '<< generating ' $@ ' >>'
	@$(SOMAKER) $(SOFLAGS) -o $@ $^ $(LIBS)
	
$(LIB_MANAGER): HKManager.o 
	@echo '<< generating ' $@ ' >>'
	@$(SOMAKER) $(SOFLAGS) -o $@ $^ $(LIBS)
		
#$(LIB_DATAMODEL): DataModelRootDict.cpp $(DATAMODEL_OBJECTS) 
$(LIB_DATAMODEL): $(DATAMODEL_OBJECTS) 
	@echo '<< generating ' $@ ' >>'
	@$(SOMAKER) $(SOFLAGS) -o $@  $(DATAMODEL_OBJECTS)  $(LIBS)
	
# DataModel objects
%.o: DataModel-lite/%.cpp DataModel-lite/%.hpp
	@echo '<< compiling DataModel object: ' $@ ' >>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<

# Generate ROOT dictionary
DataModelRootDict.cpp: $(DATAMODEL_HEADERS) DataModel-lite/LinkDef.h
	@echo '<< generating ROOT dictionary >>'
	@rootcling -f $@ -rml $(LIB_DATAMODEL) -rmf $(PCM_LIB_DATAMODEL) \
		$(INCFLAGS) $(DATAMODEL_HEADERS) DataModel-lite/LinkDef.h

DataModelRootDict.o: DataModelRootDict.cpp
	@echo '<< compiling ROOT dictionary >>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<

$(PCM_DATAMODEL): DataModelRootDict.o
	@echo '<< installing ROOT dictionary >>'
	@mv DataModelRootDict*.pcm ../lib/.

# default rules
.cc.o: 
	@echo '<< compiling' $< '>>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $<

.cpp.o: 
	@echo '<< compiling' $< '>>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c $<

%.o: %.cc %.hh
	@echo '<< compiling' $< '>>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<

%.o: %.cpp %.hpp
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
	@rm -f ./../lib/*.pcm
	@rm -f ./*Dict*
	@rm -f ./*.pcm
	
