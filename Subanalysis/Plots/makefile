CXX=g++
CXXFLAGS=`root-config --cflags --glibs`

OBJPATH=obj
SRCPATH=src
INCPATH=/internal/big_one/4/users/gamrat/scripts/Include/
BINPATH=bin

LIBS=-L$(INCPATH) -lrec -lm -llapack -lcurl

SRC=$(SRCPATH)/plots_norm.cpp plots_main.cpp
OBJ=$(OBJPATH)/plots_norm.o $(OBJPATH)/plots_main.o

PROJECT=plots.exe
RM=rm
RMFLAGS=-f

all: $(OBJ)

# $(BINPATH)/$(PROJECT): $(OBJPATH)/main.o $(OBJ)
# 	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

$(OBJPATH)/plots_norm.o: $(SRCPATH)/plots_norm.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/plots_main.o: plots_main.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJ) $(BINPATH)/$(PROJECT)