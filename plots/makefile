CXX=g++
CXXFLAGS=`root-config --cflags --glibs`
LIBS=-L/internal/big_one/4/users/gamrat/scripts/Include/ -lrec -lm -llapack
PROJECT_GEN=aux_plots
RM=rm
RMFLAGS=-f

all: $(PROJECT_GEN)

$(PROJECT_GEN): aux_plots.o
	$(CXX) $(CXXFLAGS) $(LIBS) -o $@ $^

$(PROJECT_GEN).o: aux_plots.cpp
	$(CXX) $(CXXFLAGS) $(LIBS) -c $< -o $@

clean:
	$(RM) $(RMFLAGS) *.o $(PROJECT_GEN)