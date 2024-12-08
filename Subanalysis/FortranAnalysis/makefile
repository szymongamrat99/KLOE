CFF = gfortran
LDFLAGS = -lgfortran -llapack -lblas -lkernlib 
CXXFLAGS = -fPIC `root-config --cflags --glibs` -lm -MD -MP

VPATH = src:obj
OBJDIR = obj
MODDIR = mod
SRCDIR = src
BINDIR = bin

LIBS=-lcurl -lpdg -lfort

all: $(BINDIR)/klspm00_analysis.exe $(OBJDIR)/test_api.o $(OBJDIR)/analysis_library.o $(OBJDIR)/analysismodule.o $(MODDIR)/analysismodule.mod

$(BINDIR)/klspm00_analysis.exe: $(OBJDIR)/test_api.o $(OBJDIR)/analysis_library.o
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(OBJDIR)/test_api.o: $(SRCDIR)/test_api.cpp
	$(CXX) $(LIBS) -c $^ -o $@

$(OBJDIR)/analysismodule.o: $(SRCDIR)/analysismodule.f
	$(CFF) -c $< -o $@ -J$(MODDIR)

$(OBJDIR)/analysis_library.o: $(SRCDIR)/analysis_library.f $(OBJDIR)/analysismodule.o
	$(CFF) -c $< -o $@ -J$(MODDIR)

clean:
	rm -f $(BINDIR)/klspm00_analysis.exe $(OBJDIR)/*.o $(MODDIR)/*.mod