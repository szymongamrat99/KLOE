# STANDARD DIRECTORIES FOR MAIN COMPILATION
MAKE_DIR = $(PWD)
SUBANALYSIS_DIR = $(MAKE_DIR)/Subanalysis

# DIRECTORIES FOR SUBCOMPILATIONS
CPFIT_DIR 		 := $(SUBANALYSIS_DIR)/CPFit
EFFICIENCY_DIR := $(SUBANALYSIS_DIR)/EfficiencyAnalysis
FORTRAN_DIR    := $(SUBANALYSIS_DIR)/FortranAnalysis
GENERATED_DIR  := $(SUBANALYSIS_DIR)/GeneratedVars
OMEGAREC_DIR   := $(SUBANALYSIS_DIR)/OmegaRec 
NEUTREC_DIR    := $(SUBANALYSIS_DIR)/Neutrec

# WORK DIRECTORIES
OBJDIR := obj
BINDIR := bin

OBJ :=
OBJ += $(OBJDIR)/CPVAnalysis.o

MAKE_CPFIT := $(CPFIT_DIR)/$(OBJDIR)/cpfit_main.o
MAKE_CPFIT += $(CPFIT_DIR)/$(OBJDIR)/cp_fit_mc_data.o

LIBS= -lm -llapack -lrec -lfort -lcurl

PROJECT := $(BINDIR)/KLSPM00.exe

CXX = g++
CXXFLAGS := `root-config --cflags --glibs` -g3

all: $(OBJ) $(MAKE_CPFIT) $(PROJECT)

$(PROJECT): $(OBJ) $(MAKE_CPFIT)
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@
	
$(OBJDIR)/CPVAnalysis.o: CPVAnalysis.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(MAKE_CPFIT):
	@$(MAKE) -C $(CPFIT_DIR)

clean: 
	@$(MAKE) -C $(CPFIT_DIR) clean
	rm -f $(OBJDIR)/*.o





