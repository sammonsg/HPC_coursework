#include <cstdlib>
#include <iostream>

using namespace std;

int addr(int r, int c, int cols){
    // Helper function to index from human-indices to c++ vector notation
    return (r-1)*cols+(c-1);
}

int b_addr(int r, int c, int cols){
    // Helper function to index from human-indices to c++ vector notation
    int col = 1 + (cols+1)/2 + c - r;
    cout << "Original coord: " << r << "-" << c;
    cout << "   Banded column: " << col << endl;
    return (r-1)*cols+(c-1);
}

int addr(int r){
    // Helper function to index from human-indices to c++ vector notation
    return r-1;
}

void print_m(double* M, int rows){
    for (int c = 1; c < rows + 1; c++){
        for (int r = 1; r < rows + 1; r++){
            cout << M[addr(c, r, rows)] << "\t";
        }
        cout << endl << endl;
    }
}

void print_v(double* M, int rows){
    for (int r = 1; r < rows + 1; r++){
        cout << M[addr(r)] << endl;
    }
}
