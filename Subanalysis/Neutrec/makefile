CXX=g++
CXXFLAGS=`root-config --cflags --glibs`

OBJPATH=obj
SRCPATH=src
INCPATH=/internal/big_one/4/users/gamrat/scripts/Include/
BINPATH=bin

LIBS=-L$(INCPATH) -lrec -lm -llapack -lcurl

SRC=$(SRCPATH)/trilateration_neu_rec.cpp $(SRCPATH)/trilateration_neu_rec_kin_fit_aux.cpp $(SRCPATH)/triangle_neu_rec.cpp $(SRCPATH)/comp_of_methods.cpp
OBJ=$(OBJPATH)/trilateration_neu_rec.o $(OBJPATH)/trilateration_neu_rec_kin_fit_aux.o $(OBJPATH)/triangle_neu_rec.o $(OBJPATH)/comp_of_methods.o

PROJECT=tri_rec.exe

RM=rm
RMFLAGS=-f

all: $(OBJ) $(BINPATH)/$(PROJECT)

$(BINPATH)/$(PROJECT): $(OBJPATH)/main.o $(OBJ)
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

$(OBJPATH)/trilateration_neu_rec.o : $(SRCPATH)/trilateration_neu_rec.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/triangle_neu_rec.o : $(SRCPATH)/triangle_neu_rec.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/trilateration_neu_rec_kin_fit_aux.o: $(SRCPATH)/trilateration_neu_rec_kin_fit_aux.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/comp_of_methods.o : $(SRCPATH)/comp_of_methods.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/main.o: neutrec_main.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJ) $(BINPATH)/$(PROJECT)