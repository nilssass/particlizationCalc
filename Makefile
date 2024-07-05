ODIR           = obj

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXX           = /opt/homebrew/Cellar/llvm/17.0.4/bin/clang++
CXXFLAGS      = -fPIC -O3 -fopenmp  # Add -fopenmp for OpenMP support
LD            = /opt/homebrew/Cellar/llvm/17.0.4/bin/clang++
LDFLAGS       = -O3 -fopenmp -lomp # Add -fopenmp for OpenMP support
FFLAGS        = -fPIC $(ROOTCFLAGS) -O3

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

_HYDROO        = DecayChannel.o ParticlePDG2.o DatabasePDG2.o UKUtility.o gen.o \
                particle.o main.o interpolation.o
 
# VPATH = src:../UKW
HYDROO = $(patsubst %,$(ODIR)/%,$(_HYDROO))

TARGET = calc
#------------------------------------------------------------------------------

$(TARGET): $(HYDROO)
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"

clean:
		@rm -f $(ODIR)/*.o $(TARGET)

$(ODIR)/%.o: src/%.cpp src/const.h
		$(CXX) $(CXXFLAGS) -c $< -o $@

