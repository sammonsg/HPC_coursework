#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

int b_addr(int r, int c, int bw){
    // Helper function to index from human-indices to c++ vector notation
    int rows = 3 * bw + 1;
    int row = 1 + 2 * bw - c + r;
    int ret = (c-1)*rows+(row-1);
    return ret;
}

void print_banded_m(double* M, int rows, int cols){
    // Transposing is necessary to print the matrix properly on screen
    cout << endl << cols << " by " << rows << " matrix, transposed:" << endl;
    for (int r = 1; r < rows + 1; r++){
        for (int c = 1; c < cols + 1; c++){
            cout << setw(13) << M[(c-1)*rows+(r-1)];
        }
        cout << endl << endl;
    }
}

void clone_vector(double* master, double* target, int length){
    for (int a = 0; a < length; a++){
        target[a] = master[a];
    }
}

void shift_vec(double*& pres, double*& futu){
    double* temp = pres;
    pres = futu;
    futu = temp;
    return;
}

void shift_vec(double*& past, double*& pres, double*& futu){
    double* temp = past;
    past = pres;
    pres = futu;
    futu = temp;
    return;
}

void print_v(double* M, int cols){
    cout << endl << "1 by " << cols << " vector, transposed:" << endl;
    for (int c = 0; c < cols; c++){
        cout << setw(14) << M[c] << endl;
    }
}

void print_pos_v(double* M, int cols){
    if (cols % 3 != 0){
        cout << "Cannot print" << endl;
        return;
    }
    int rows = cols / 3;
    cout << endl << "1 by " << cols << " vector, in coordinate form" << endl;
    cout << setw(16) << "Node" << setw(16) << "x-disp" << setw(16) << "y-disp" << setw(16) << "noderotation" << endl;
    for (int r = 0; r < rows; r++){
        cout << setw(16) << r << setw(16) << M[3*r] << setw(16) << M[3*r+1] << setw(16) << M[3*r+2] << endl;
    }
}

void m_diag_add(double* L, double* R, double r_fact, int l, int bw){
    int rows = 3 * bw + 1;
    for (int a = 0; a < l; a++){
        L[a * rows + 2 * bw] += R[a] * r_fact;
    }
}
