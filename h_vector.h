//Joe Snider
//12/05
//
//Derive a vector class that can be instantiated with an int.
//Just creates a blank vector when an int is passed to the constructor.
//This allows compatibility with the histogram class.

#include <vector>

using namespace std;

template<class T>
class CHVector : public vector<T>{
public:
   CHVector(const int& inInt):vector<T>() {}
};
