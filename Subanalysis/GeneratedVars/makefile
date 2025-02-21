CXX=g++
CXXFLAGS=`root-config --cflags --glibs`

OBJPATH=obj
SRCPATH=src
INCPATH=/internal/big_one/4/users/gamrat/scripts/Include/
BINPATH=bin

LIBS=-L$(INCPATH) -lrec -lm -llapack -lcurl

SRC=$(SRCPATH)/split_channels.cpp genvars_main.cpp #$(SRCPATH)/generated_variables.cpp
OBJ=$(OBJPATH)/split_channels.o $(OBJPATH)/genvars_main.o #$(OBJPATH)/generated_variables.o 

PROJECT=gen_vars.exe
RM=rm
RMFLAGS=-f

all: $(OBJ)

# $(BINPATH)/$(PROJECT): $(OBJPATH)/genvars_main.o $(OBJ)
# 	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

$(OBJPATH)/split_channels.o: $(SRCPATH)/split_channels.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@	

# $(OBJPATH)/generated_variables.o: $(SRCPATH)/generated_variables.cpp
# 	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/genvars_main.o: genvars_main.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJ) $(BINPATH)/$(PROJECT)