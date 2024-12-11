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
SRCDIR := src

OBJ :=
OBJ += $(OBJDIR)/CPVAnalysis.o

HEADER := ${WORKDIR}/scripts/Include/const.h

MAKE_CPFIT := $(CPFIT_DIR)/$(OBJDIR)/cpfit_main.o
MAKE_CPFIT += $(CPFIT_DIR)/$(OBJDIR)/cp_fit_mc_data.o

SRC_CPFIT := $(CPFIT_DIR)/$(SRCDIR)/cp_fit_mc_data.cpp
SRC_CPFIT += $(CPFIT_DIR)/cpfit_main.cpp

MAKE_GENVARS := $(GENERATED_DIR)/$(OBJDIR)/genvars_main.o
MAKE_GENVARS += $(GENERATED_DIR)/$(OBJDIR)/generated_variables.o
MAKE_GENVARS += $(GENERATED_DIR)/$(OBJDIR)/split_channels.o

SRC_GENVARS := $(GENERATED_DIR)/$(SRCDIR)/split_channels.cpp
SRC_GENVARS += $(GENERATED_DIR)/$(SRCDIR)/generated_variables.cpp
SRC_GENVARS += $(GENERATED_DIR)/genvars_main.cpp

MAKE_OMEGAREC := $(OMEGAREC_DIR)/$(OBJDIR)/omegarec_main.o
MAKE_OMEGAREC += $(OMEGAREC_DIR)/$(OBJDIR)/omega_rec_kin_fit.o
MAKE_OMEGAREC += $(OMEGAREC_DIR)/$(OBJDIR)/plots.o

SRC_OMEGAREC := $(OMEGAREC_DIR)/$(SRCDIR)/omega_rec_kin_fit.cpp
SRC_OMEGAREC += $(OMEGAREC_DIR)/$(SRCDIR)/plots.cpp
SRC_OMEGAREC += $(OMEGAREC_DIR)/omegarec_main.cpp

LIBS= -lm -llapack -lrec -lfort -lcurl

PROJECT := $(BINDIR)/KLSPM00.exe

CXX = g++
CXXFLAGS := `root-config --cflags --glibs` -g3

all: $(OBJ) $(MAKE_CPFIT) $(MAKE_GENVARS) $(MAKE_OMEGAREC) $(PROJECT) $(SRC_CPFIT) $(SRC_OMEGAREC) $(SRC_GENVARS)

$(PROJECT): $(OBJ) $(MAKE_CPFIT) $(MAKE_GENVARS) $(MAKE_OMEGAREC)
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@
	
$(OBJDIR)/CPVAnalysis.o: CPVAnalysis.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(MAKE_CPFIT):
	@$(MAKE) -C $(CPFIT_DIR)

$(MAKE_OMEGAREC):
	@$(MAKE) -C $(OMEGAREC_DIR)

$(MAKE_GENVARS):
	@$(MAKE) -C $(GENERATED_DIR)

clean: 
	@$(MAKE) -C $(CPFIT_DIR) clean
	@$(MAKE) -C $(GENERATED_DIR) clean
	@$(MAKE) -C $(OMEGAREC_DIR) clean
	rm -f $(OBJDIR)/*.o





