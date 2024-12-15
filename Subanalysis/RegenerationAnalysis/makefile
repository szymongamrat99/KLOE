CXX=g++
CXXFLAGS=`root-config --cflags --glibs`

OBJPATH=obj
SRCPATH=src
INCPATH=/internal/big_one/4/users/gamrat/scripts/Include/
BINPATH=bin

LIBS=-L$(INCPATH) -lrec -lm -llapack -lcurl

SRC=$(SRCPATH)/regen_cut_analysis.cpp $(SRCPATH)/plots.cpp regen_main.cpp
OBJ=$(OBJPATH)/regen_cut_analysis.o $(OBJPATH)/plots.o $(OBJPATH)/regen_main.o

PROJECT=regen_cut_analysis.exe
RM=rm
RMFLAGS=-f

all: $(OBJ) $(SRC)

# $(BINPATH)/$(PROJECT): $(OBJPATH)/main.o $(OBJ)
# 	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

$(OBJPATH)/plots.o: $(SRCPATH)/plots.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/regen_cut_analysis.o: $(SRCPATH)/regen_cut_analysis.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/regen_main.o: regen_main.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJ) $(BINPATH)/$(PROJECT)