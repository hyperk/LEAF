
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

#all: .depend $(DICT).cc libBQFitter.so $(OBJECT)
all: .depend ../lib/libBQFitter.so

# library
../lib/libBQFitter.so: BQFitter.o 

#$(DICT).o
	@echo '<< generating ../lib/libBQFitter.so >>'
#	@$(SOMAKER) $(SOFLAGS) -o ../lib/libBQFitter.so  BQFitter.o $(DICT).o $(LIBS) 
	@$(SOMAKER) $(SOFLAGS) -o ../lib/libBQFitter.so  BQFitter.o $(LIBS) 
	
# root class dictionary
#$(DICT).cc:  LinkDef.h
#	@echo '<< generating CINT >>'
#	@rootcint -f $(DICT).cc -c $(INCFLAGS) $(HEADERS) LinkDef.h

	
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
	
