CXX=g++
CXXFLAGS=`root-config --cflags --glibs`

OBJPATH=obj
SRCPATH=src
INCPATH=/internal/big_one/4/users/gamrat/scripts/Include/
BINPATH=bin

LIBS=-L$(INCPATH) -lrec -lm -llapack -lcurl

SRC=$(SRCPATH)/kchrec_Kmass.cpp KchRec_main.cpp
OBJ=$(OBJPATH)/kchrec_Kmass.o $(OBJPATH)/KchRec_main.o

PROJECT=Kchrec.exe

RM=rm
RMFLAGS=-f

all: $(OBJ)

# $(BINPATH)/$(PROJECT): $(OBJPATH)/Kchrec_main.o $(OBJ)
# 	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

$(OBJPATH)/kchrec_Kmass.o : $(SRCPATH)/kchrec_Kmass.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/KchRec_main.o: KchRec_main.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJ) $(BINPATH)/$(PROJECT)