#Env vars
CXX = g++
CXXFLAGS = -fPIC `root-config --cflags --glibs` -lm
LDFLAGS = -shared -g

kloe_library.so: kloe_class.o
	$(CXX) -o kloe_library.so $(LDFLAGS) kloe_class.o

kloe_class.o: kloe_class.cpp kloe_class.h
	$(CXX) $(CXXFLAGS) -c kloe_class.cpp -o kloe_class.o  
