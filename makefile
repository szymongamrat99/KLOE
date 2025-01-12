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
REGEN_DIR    	 := $(SUBANALYSIS_DIR)/RegenerationAnalysis
PLOTS_DIR    	 := $(SUBANALYSIS_DIR)/Plots

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
MAKE_OMEGAREC += $(OMEGAREC_DIR)/$(OBJDIR)/omega_rec.o
MAKE_OMEGAREC += $(OMEGAREC_DIR)/$(OBJDIR)/plots_df.o

SRC_OMEGAREC := $(OMEGAREC_DIR)/$(SRCDIR)/omega_rec.cpp
SRC_OMEGAREC += $(OMEGAREC_DIR)/$(SRCDIR)/plots_df.cpp
SRC_OMEGAREC += $(OMEGAREC_DIR)/omegarec_main.cpp

MAKE_REGEN := $(REGEN_DIR)/$(OBJDIR)/regen_main.o
MAKE_REGEN += $(REGEN_DIR)/$(OBJDIR)/regen_cut_analysis.o
MAKE_REGEN += $(REGEN_DIR)/$(OBJDIR)/plots.o

SRC_REGEN := $(REGEN_DIR)/$(SRCDIR)/regen_cut_analysis.cpp
SRC_REGEN := $(REGEN_DIR)/$(SRCDIR)/plots.cpp
SRC_REGEN += $(REGEN_DIR)/regen_main.cpp

MAKE_PLOTS := $(PLOTS_DIR)/$(OBJDIR)/plots_main.o
MAKE_PLOTS += $(PLOTS_DIR)/$(OBJDIR)/plots_norm.o

SRC_PLOTS := $(PLOTS_DIR)/$(SRCDIR)/plots_norm.cpp
SRC_PLOTS += $(PLOTS_DIR)/plots_main.cpp

LIBS= -lm -llapack -lrec -lfort -lcurl

PROJECT := $(BINDIR)/KLSPM00.exe

CXX = g++
CXXFLAGS := `root-config --cflags --glibs` -g3 -fopenmp

all: $(OBJ) $(MAKE_CPFIT) $(MAKE_GENVARS) $(MAKE_OMEGAREC) $(MAKE_REGEN) $(MAKE_PLOTS) $(PROJECT)

$(PROJECT): $(OBJ) $(MAKE_CPFIT) $(MAKE_GENVARS) $(MAKE_OMEGAREC) $(MAKE_REGEN) $(MAKE_PLOTS) $(HEADER)
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@
	
$(OBJDIR)/CPVAnalysis.o: CPVAnalysis.cpp $(HEADER)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(MAKE_CPFIT): $(SRC_CPFIT) $(HEADER)
	@$(MAKE) -C $(CPFIT_DIR)

$(MAKE_OMEGAREC): $(SRC_OMEGAREC) $(HEADER)
	@$(MAKE) -C $(OMEGAREC_DIR)

$(MAKE_GENVARS): $(SRC_GENVARS) $(HEADER)
	@$(MAKE) -C $(GENERATED_DIR)

$(MAKE_REGEN): $(SRC_REGEN) $(HEADER)
	@$(MAKE) -C $(REGEN_DIR)

$(MAKE_PLOTS): $(SRC_PLOTS) $(HEADER)
	@$(MAKE) -C $(PLOTS_DIR)

clean: 
	@$(MAKE) -C $(CPFIT_DIR) clean
	@$(MAKE) -C $(GENERATED_DIR) clean
	@$(MAKE) -C $(OMEGAREC_DIR) clean
	@$(MAKE) -C $(REGEN_DIR) clean
	rm -f $(OBJDIR)/*.o





