
OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)
DICT  		= BenjaminFitterCint

#SOURCES 	= $(wildcard *.cc)
#HEADERS 	= $(wildcard *.hh)

SOURCES 	= BQFitter.cc
HEADERS 	= BQFitter.hh

include ../Makefile/Makefile.${OSNAME}

# set compiler options for ROOT
CXXFLAGS += $(shell root-config --cflags)
CXXFLAGS += '-fPIC' -std=c++11

INCFLAGS = -I.
INCFLAGS += -I$(WCSIMDIR)/include

LIBS += ${WCSIMDIR}/libWCSimRoot.so 
LIBS += $(shell root-config --libs) -lMinuit 

LIB_MANAGER = ../lib/libHKManager.so
LIB_ANALYSIS = ../lib/libHKAstroAnalysis.so
LIB_FITTER = ../lib/libBQFitter.so

LIBS_OBJECTS = ../lib/libHKManager.so ../lib/libHKAstroAnalysis.so ../lib/libBQFitter.so

#all: .depend $(DICT).cc libBQFitter.so $(OBJECT)
all: .depend $(LIBS_OBJECTS)

# library
../lib/libBQFitter.so: BQFitter.o 
	@if [ ! -d ../lib ]; then \
		mkdir ../lib; \
	fi
	@echo '<< generating ../lib/libBQFitter.so >>'
	@$(SOMAKER) $(SOFLAGS) -o ../lib/libBQFitter.so  BQFitter.o $(LIBS) $(LIB_MANAGER) $(LIB_ANALYSIS)
	
../lib/libHKManager.so: HKManager.o 
	@echo '<< generating ../lib/libHKManager.so >>'
	@if [ ! -d ../lib ]; then \
		mkdir ../lib; \
	fi
	@$(SOMAKER) $(SOFLAGS) -o ../lib/libHKManager.so  HKManager.o $(LIBS) 
	
../lib/libHKAstroAnalysis.so: HKAstroAnalysis.o 
	@echo '<< generating ../lib/libHKManager.so >>'
	@if [ ! -d ../lib ]; then \
		mkdir ../lib; \
	fi
	@$(SOMAKER) $(SOFLAGS) -o ../lib/libHKAstroAnalysis.so  HKAstroAnalysis.o $(LIBS) $(LIB_MANAGER)

	
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
	@rm -f ./$(DICT)*
	
