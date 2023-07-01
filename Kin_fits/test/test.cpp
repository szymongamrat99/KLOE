#include <iostream>

extern "C" void print_hi_(void);

using namespace std;

int main() 
{
  print_hi_();
  cout << "Hello from C++" << endl;
  return 0;
}