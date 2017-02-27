#include <iostream>
#include <cstdlib>

using namespace std;

int main(int argc, char* argv[]) {
    cout << "argc = " << argc << endl;
    for(int i = 1; i < argc; i++)
        cout << atoi(argv[i])*2 << endl;
    return 0;
}