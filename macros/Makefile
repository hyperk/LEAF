CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated -fpermissive -std=c++11

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

LIBS 		= $(ROOTGLIBS) -lMinuit $(WCSIMDIR)/libWCSimRoot.so

INC = $(WCSIMDIR)/include
SRC= $(WCSIMDIR)/src

CXXFLAGS += -I$(LIBS) -I$(SRC) -I$(INC) 

TARGET= LEAFOutputAnalysisHybrid_leafclass LEAFOutputAnalysisHybrid_master AnalyzeWSHierarchy ProducePDF #FitVertexLE

all: $(TARGET)
LEAFOutputAnalysisHybrid_leafclass: LEAFOutputAnalysisHybrid_leafclass.o
LEAFOutputAnalysisHybrid_master: LEAFOutputAnalysisHybrid_master.o
ProducePDF: ProducePDF.o
AnalyzeWSHierarchy: AnalyzeWSHierarchy.o
FitVertexLE: FitVertexLE.o

%: %.o
	@echo "Now make $@"
	@$(CPP) -o $@ $< $(CXXFLAGS) $(LIBS) 
	@echo "..Compile done! "

%.o: %.c
	@echo "$<"
	@echo "Start Compiling $<"
	@$(CPP) $(CXXFLAGS) -c $<
	@echo ".. Compiling Object Files $<   --> done"
	@echo ""

%.o: %.cpp
	@echo "$<"
	@echo "Start Compiling $<"
	@$(CPP) $(CXXFLAGS) -c $<
	@echo ".. Compiling Object Files $<   --> done"
	@echo ""

clean: 
	@echo "Now Clean Up"
	rm -f $(TARGET) *~ *.o *.o~ core
