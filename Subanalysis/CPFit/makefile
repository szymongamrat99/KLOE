CXX=g++
CXXFLAGS=`root-config --cflags --glibs` -g3

OBJPATH=obj
SRCPATH=src
BINPATH=bin
INCPATH=inc
CONSTPATH=/internal/big_one/4/users/gamrat/scripts/Include

LIBS= -lm -llapack -lrec -lcurl

SRC=$(SRCPATH)/cp_fit_mc_data.cpp cpfit_main.cpp
OBJ=$(OBJPATH)/cp_fit_mc_data.o $(OBJPATH)/cpfit_main.o
INC=$(INCPATH)/cpfit.hpp $(CONSTPATH)/const.h

PROJECT=cpv_fit.exe

RM=rm
RMFLAGS=-f

all: $(OBJ) $(BINPATH)/$(PROJECT)

#$(BINPATH)/$(PROJECT): $(OBJPATH)/main.o $(OBJ)
#	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

$(OBJPATH)/cp_fit_mc_data.o: $(SRCPATH)/cp_fit_mc_data.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/cpfit_main.o: cpfit_main.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJ)