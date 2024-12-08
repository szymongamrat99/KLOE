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

OBJ :=
OBJ += $(OBJDIR)/CPVAnalysis.o

COLOR_ON = color

CXX = g++
CXXFLAGS := `root-config --cflags --glibs` -g3

#export MAKE_DIR CXX CXXFLAGS

$(OBJPATH)/CPVAnalysis.o: CPVAnalysis.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

#	@$(MAKE) -C $(CPFIT_DIR)

.PHONY: 
	clean

clean: 
	@$(MAKE) -C $(CPFIT_DIR) clean





