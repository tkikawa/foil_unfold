CPP 		= g++
CXXFLAGS	= -g -O3 -Wall -fPIC -D_REENTRANT -Wno-deprecated

ROOTCFLAGS	:= $(shell root-config --cflags)
ROOTLIBS     	:= $(shell root-config --libs)
ROOTGLIBS    	:= $(shell root-config --glibs)
CXXFLAGS	+= $(ROOTCFLAGS)

LIBS 		= $(ROOTLIBS)

CPPSRC = Fitter.cc Material.cc Foil.cc Cover.cc
CPPOBJ = $(CPPSRC:.cc=.o)

TARGET=Unfold
OBJS=Unfold.o

$(TARGET): $(OBJS) $(CPPOBJ)
	@echo "Now making $@"
	@$(CPP) -o $@ $(OBJS) $(CPPOBJ) $(CXXFLAGS) $(LIBS)
	@echo "Compile done! \(^o^)/"

.cc.o:
	@echo "Compiling object file $<"
	@$(CPP) -c $(CXXFLAGS) $<

clean: 
	@echo "Now cleaning up"
	rm -f $(TARGET) *~ *.o *.o~ core
