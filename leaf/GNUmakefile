
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

SOURCES 	= LEAF.cc
HEADERS 	= LEAF.hh

include ../Makefile/Makefile.${OSNAME}


# set compiler options for ROOT
CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += '-fPIC' -std=c++11

INCFLAGS = -I.
INCFLAGS += -I$(WCSIMDIR)/include
INCFLAGS += -I./DataModel

LIBS += ${WCSIMDIR}/libWCSimRoot.so 
LIBS += $(shell root-config --libs) -lMinuit 

LIB_DATAMODEL = ../lib/libDataModel.so
LIB_MANAGER = ../lib/libHKManager.so
LIB_FITTER = ../lib/libLEAF.so
LIB_ASTROANALYSIS = ../lib/libHKAstroAnalysis.so

LIBS_OBJECTS = $(LIB_DATAMODEL) $(LIB_MANAGER) $(LIB_FITTER) #$(LIB_ASTROANALYSIS)
DATAMODEL_OBJECTS = GeoTools.o TimeDelta.o HitInfo.o PMTInfo.o Geometry.o 


all: .depend $(LIBS_OBJECTS)



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
	
$(LIB_DATAMODEL): $(DATAMODEL_OBJECTS)
	@echo '<< generating ' $@ ' >>'
	@$(SOMAKER) $(SOFLAGS) -o $@  $(DATAMODEL_OBJECTS)  $(LIBS)
	

# DataModel objects
GeoTools.o: DataModel/GeoTools.cpp DataModel/GeoTools.h
	@echo '<< compiling DataModel object: ' $@ ' >>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<
	
TimeDelta.o: DataModel/TimeDelta.cpp DataModel/TimeDelta.h
	@echo '<< compiling DataModel object: ' $@ ' >>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<
	
HitInfo.o: DataModel/HitInfo.cpp DataModel/HitInfo.h
	@echo '<< compiling DataModel object: ' $@ ' >>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<
	
PMTInfo.o: DataModel/PMTInfo.cpp DataModel/PMTInfo.h
	@echo '<< compiling DataModel object: ' $@ ' >>'
	@$(CXX) $(CXXFLAGS) $(INCFLAGS) -c -o $@ $<
	
Geometry.o: DataModel/Geometry.cpp DataModel/Geometry.h
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
	@rm -f ./*.o
	@rm -f ./../lib*.so
	
