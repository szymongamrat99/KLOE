#include <iostream>


using namespace std;

int main()
{

  try
  {
    int age = 0;

    cout << "Input Your age: ";
    cin >> age;
    cout << endl;
    if(age >= 18)
    {
      cout << "Access granted!" << endl;
    }
    
  }
  catch(float myNum)
  {
    cerr << "You are not old enough!" << endl;
    cerr << "Age is: " << myNum;
  }
  catch(...)
  {
    cerr << "unknown exception" << endl;
  }
  

  return 0;
}