CXX=g++
CXXFLAGS=`root-config --cflags --glibs`

OBJPATH=obj
SRCPATH=src
INCPATH=/internal/big_one/4/users/gamrat/scripts/Include/
BINPATH=bin

LIBS=-L$(INCPATH) -lrec -lm -llapack -lcurl

SRC=$(SRCPATH)/trilateration_neu_rec_kin_fit_aux.cpp neutrec_main.cpp $(SRCPATH)/triangle_neu_rec.cpp $(SRCPATH)/comp_of_methods.cpp #$(SRCPATH)/trilateration_neu_rec.cpp $(SRCPATH)/trilateration_neu_rec_kin_fit_aux.cpp neutrec_main.cpp
OBJ=$(OBJPATH)/trilateration_neu_rec_kin_fit_aux.o $(OBJPATH)/neutrec_main.o $(OBJPATH)/triangle_neu_rec.o $(OBJPATH)/comp_of_methods.o #$(OBJPATH)/trilateration_neu_rec.o $(OBJPATH)/trilateration_neu_rec_kin_fit_aux.o $(OBJPATH)/neutrec_main.o

PROJECT=tri_rec.exe

RM=rm
RMFLAGS=-f

all: $(OBJ)

# $(BINPATH)/$(PROJECT): $(OBJPATH)/neutrec_main.o $(OBJ)
# 	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

# $(OBJPATH)/trilateration_neu_rec.o : $(SRCPATH)/trilateration_neu_rec.cpp
# 	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/triangle_neu_rec.o : $(SRCPATH)/triangle_neu_rec.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/trilateration_neu_rec_kin_fit_aux.o: $(SRCPATH)/trilateration_neu_rec_kin_fit_aux.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/comp_of_methods.o : $(SRCPATH)/comp_of_methods.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/neutrec_main.o: neutrec_main.cpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJ) $(BINPATH)/$(PROJECT)