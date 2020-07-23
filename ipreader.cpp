
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main () {
  string line;
  ifstream iplist ("./config/iplist.ini");
  if (iplist.is_open())
  {
      while ( getline (iplist,line) )
      {
        cout << line << '\n';
      }
      iplist.close();
  }

  else cout << "Unable to open file"; 

  return 0;
}