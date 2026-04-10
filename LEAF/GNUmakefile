
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
INCFLAGS += -I../DataManager/DataModel-lite

LIBS += -L../libWCSIM -lWCSimRoot
LIBS += $(shell root-config --libs) -lMinuit 

LIB_FITTER = ../lib/libLEAF.so

LIBS_OBJECTS = $(LIB_FITTER) 

all: .depend ROOT_CHECK DIRECTORY_CHECK $(LIBS_OBJECTS) 

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
$(LIB_FITTER): LEAF.o LEAFConfig.o LEAFUtilities.o SplineUtilities.o
	@echo '<< generating ' $@ ' >>'
	@$(SOMAKER) $(SOFLAGS) -o $@ $^ $(LIBS)
	
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
	@rm -f $(LIB_FITTER)
	
