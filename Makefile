ODIR           = obj

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXX           = /opt/homebrew/Cellar/llvm/18.1.5/bin/clang++
# CXXFLAGS      = -g -fPIC -O3 -fopenmp -Xpreprocessor # Add -fopenmp for OpenMP support
CXXFLAGS      = -fPIC -O3 -fopenmp -Xpreprocessor # Add -fopenmp for OpenMP support
LD            = /opt/homebrew/Cellar/llvm/18.1.5/bin/clang++
LDFLAGS       = -O3 -fopenmp -lomp # Add -fopenmp for OpenMP support
# FFLAGS        = -fPIC $(ROOTCFLAGS) -O3
FFLAGS        = -fPIC -O3

CXXFLAGS     += $(ROOTCFLAGS)
# LIBS          = $(ROOTLIBS) $(SYSLIBS)
# GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

LIBS          = $(SYSLIBS)
GLIBS         = $(SYSLIBS)

# _HYDROO        = DecayChannel.o ParticlePDG2.o DatabasePDG2.o UKUtility.o gen.o \
#                 particle.o main.o interpolation.o utils.o elements.o
_HYDROO        =  main.o utils.o element.o engine.o pdg_particle.o surface.o
 
# VPATH = src:../UKW
HYDROO = $(patsubst %,$(ODIR)/%,$(_HYDROO))

TARGET = calc
#------------------------------------------------------------------------------

$(TARGET): $(HYDROO)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"

clean:
		@rm -f $(ODIR)/*.o $(TARGET)

$(ODIR)/%.o: src/%.cpp 
		$(CXX) $(CXXFLAGS) -c $< -o $@ 

$(ODIR)/%.o: src/%.hpp 
		$(CXX) $(CXXFLAGS) -c $< -o $@ 

# $(ODIR)/cmdparser.h.pch: src/cmdparser.hpp 
# 		$(CXX) -cc1 cmdparser.hpp -emit-pch -o cmdparser.h.pch
