CXX=g++
CXXFLAGS=`root-config --cflags --glibs`

OBJPATH=obj
SRCPATH=src
INCPATH=/internal/big_one/4/users/gamrat/scripts/Include/
BINPATH=bin

LIBS=-L$(INCPATH) -lrec -lm -llapack

SRC=$(SRCPATH)/cp_fit_mc_data.cpp
OBJ=$(OBJPATH)/cp_fit_mc_data.o

PROJECT=cpv_fit.exe

RM=rm
RMFLAGS=-f

all: $(OBJ) $(BINPATH)/$(PROJECT)

$(BINPATH)/$(PROJECT): $(OBJPATH)/main.o $(OBJ)
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

$(OBJPATH)/cp_fit_mc_data.o: $(SRCPATH)/cp_fit_mc_data.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/main.o: cpfit_main.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJ) $(BINPATH)/$(PROJECT)