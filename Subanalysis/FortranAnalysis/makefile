CFF = gfortran
LDFLAGS = -lgfortran -llapack -lblas -lkernlib -Wall
CXXFLAGS = `root-config --cflags --glibs` -lm -MD -MP -Wall

VPATH = src:obj
OBJDIR = obj
FORTDIR = /internal/big_one/4/users/gamrat/scripts/Include/FortranAnalysis/obj
FORTMOD = /internal/big_one/4/users/gamrat/scripts/Include/FortranAnalysis/mod
FORTSRC = /internal/big_one/4/users/gamrat/scripts/Include/FortranAnalysis/src
MODDIR = mod
SRCDIR = src
BINDIR = bin

LIBS= -L/internal/big_one/4/users/gamrat/scripts/Include/FortranAnalysis/lib/libfort.a -lcurl -lpdg

all: $(BINDIR)/klspm00_analysis.exe $(OBJDIR)/charged_tracks.o $(FORTMOD)/analysismodule.mod

$(BINDIR)/klspm00_analysis.exe: $(OBJDIR)/charged_tracks.o $(FORTDIR)/analysis_library.o 
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(OBJDIR)/charged_tracks.o: $(SRCDIR)/charged_tracks.cpp
	$(CXX) $(LIBS) -c $^ -o $@

$(FORTDIR)/analysismodule.o: $(FORTSRC)/analysismodule.f
	$(CFF) -c $< -o $@ -J$(MODDIR)

$(FORTDIR)/analysis_library.o: $(FORTSRC)/analysis_library.f $(FORTDIR)/analysismodule.o
	$(CFF) -c $< -o $@ -J$(MODDIR)

clean:
	rm -f $(BINDIR)/klspm00_analysis.exe $(OBJDIR)/*.o $(MODDIR)/*.mod