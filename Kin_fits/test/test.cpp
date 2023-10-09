#include <iostream>

extern "C"{
  void print_hi_(int M, double *P);
}

using namespace std;

int main() 
{
  int M = 6;

  double P[36];
  print_hi_(M,P);
  cout << "Hello from C++" << endl;
  return 0;
}