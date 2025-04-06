CXX=g++
CXXFLAGS=`root-config --cflags --glibs` -g3

OBJPATH=obj
SRCPATH=src
BINPATH=bin
INCPATH=inc
CONSTPATH=/internal/big_one/4/users/gamrat/scripts/Include

LIBS= -lm -llapack -lrec -lcurl -lboost_filesystem -lboost_system

SRC=$(SRCPATH)/cov_matrix_determination.cpp $(SRCPATH)/cov_matrix_determination_control_sample.cpp covmatrix_main.cpp
OBJ=$(OBJPATH)/cov_matrix_determination.o $(OBJPATH)/cov_matrix_determination_control_sample.o $(OBJPATH)/covmatrix_main.o
INC=$(INCPATH)/covmatrix.hpp $(CONSTPATH)/const.h

PROJECT=covmatrix.exe

RM=rm
RMFLAGS=-f

all: $(OBJ)

#$(BINPATH)/$(PROJECT): $(OBJPATH)/main.o $(OBJ)
#	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

$(OBJPATH)/cov_matrix_determination.o: $(SRCPATH)/cov_matrix_determination.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/cov_matrix_determination_control_sample.o: $(SRCPATH)/cov_matrix_determination_control_sample.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJPATH)/covmatrix_main.o: covmatrix_main.cpp $(INC)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) $(OBJ)