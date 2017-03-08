#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace std;

int b_addr(int r, int c, int bw){
    // Helper function to index from human-indices to c++ vector notation
    int rows = 3 * bw + 1;
    int row = 1 + 2 * bw - c + r;
    int ret = (c-1)*rows+(row-1);
    // cout << "Original coord: " << r << "-" << c;
    // cout << "   Banded row: " << row << "\tReturn value: " << ret << endl;
    return ret;
}

void print_banded_m(double* M, int rows, int cols){
    // Transposing is necessary to print the matrix properly on screen
    for (int r = 1; r < rows + 1; r++){
        for (int c = 1; c < cols + 1; c++){
            // cout << addr(r, c, rows) << endl;
            cout << setw(12) << M[(c-1)*rows+(r-1)] << "\t";
        }
        cout << endl << endl;
    }
}

void print_v(double* M, int cols){
    for (int c = 0; c < cols; c++){
        cout << M[c] << endl;
    }
}
