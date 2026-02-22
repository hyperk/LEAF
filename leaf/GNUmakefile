
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

SOURCES 	= LEAF.cc
HEADERS 	= LEAF.hh

include ../Makefile/Makefile.${OSNAME}



# set compiler options for ROOT
CXXFLAGS_BASE += $(shell root-config --cflags)
CXXFLAGS_BASE += '-fPIC' -Wall -Wpedantic -Wno-long-long -DDATAMODEL_LITE

CXXFLAGS = $(CXXFLAGS_BASE)


ROOT_VERSION = $(shell root-config --version | cut -d'.' -f1,2)
ifeq ($(shell echo "$(ROOT_VERSION) >= 6.0" | bc -l), 1)
MAKE_LEAF_WITH_ROOT6 = 1
CXXFLAGS += -std=c++17 -DROOT6 # modifies to c++17 (Nicolas)
#CXXFLAGS_ROOT += -std=c++11 -DROOT6
else
CXXFLAGS += -std=c++17 -DROOT5
endif

INCFLAGS = -I. -I$(shell root-config --incdir)
INCFLAGS += -I$(WCSIM_BUILD_DIR)/include/WCSim
INCFLAGS += -I./DataModel

LIBS += -L${WCSIM_BUILD_DIR}/lib -lWCSimRoot
LIBS += $(shell root-config --libs) -lMinuit 

LIB_DATAMODEL = ../lib/libDataModelLite.so
LIB_MANAGER = ../lib/libHKManager.so
LIB_FITTER = ../lib/libLEAF.so
LIB_ASTROANALYSIS = ../lib/libHKAstroAnalysis.so

LIBS_OBJECTS = $(LIB_DATAMODEL) $(LIB_MANAGER) $(LIB_FITTER) #$(LIB_ASTROANALYSIS)

DATAMODEL_SOURCES = $(filter-out DataModel/DataModelRootDict.cpp, $(wildcard DataModel/*.cpp))
DATAMODEL_SRC_NOP = $(DATAMODEL_SOURCES:DataModel/%=./%)

DATAMODEL_HEADERS = $(DATAMODEL_SOURCES:%.cpp=%.h)
DATAMODEL_OBJECTS = $(DATAMODEL_SRC_NOP:%.cpp=%.o) #DataModelRootDict.o


all: .depend ROOT_CHECK $(LIBS_OBJECTS)

ROOT_CHECK:
ifeq ($(MAKE_LEAF_WITH_ROOT6),1)
	@echo '<< compile with ROOT6 >>'
else
	@echo '<< compile with ROOT5 >>'
	@echo $(DATAMODEL_SOURCES)
	
endif

# library
$(LIB_FITTER): LeafConfig.o LeafInputs.o LeafUtility.o LeafSplines.o Likelihoods.o LEAF.o
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
	@if [ ! -d ../lib ]; then \
		mkdir ../lib; \
	fi
	@$(SOMAKER) $(SOFLAGS) -o $@  $(DATAMODEL_OBJECTS)  $(LIBS)
	
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
	
