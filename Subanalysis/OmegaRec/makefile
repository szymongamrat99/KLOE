CXX=g++
CXXFLAGS=`root-config --cflags --glibs`

OBJPATH=obj
SRCPATH=src
INCPATH=/internal/big_one/4/users/gamrat/scripts/Include/
BINPATH=bin

LIBS=-L$(INCPATH) -lrec -lm -llapack -lcurl

SRC=$(SRCPATH)/omega_rec.cpp $(SRCPATH)/plots.cpp omegarec_main.cpp
OBJ=$(OBJPATH)/omega_rec.o $(OBJPATH)/plots.o $(OBJPATH)/omegarec_main.o

PROJECT=omega_rec.exe
RM=rm
RMFLAGS=-f

all: $(OBJ)

# $(BINPATH)/$(PROJECT): $(OBJPATH)/main.o $(OBJ)
# 	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

$(OBJPATH)/omega_rec.o: $(SRCPATH)/omega_rec.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

#$(OBJPATH)/omega_rec_kin_fit.o: $(SRCPATH)/omega_rec_kin_fit.cpp
#	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/plots.o: $(SRCPATH)/plots.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/omegarec_main.o: omegarec_main.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJ) $(BINPATH)/$(PROJECT)