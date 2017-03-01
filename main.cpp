#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "cblas.h"

void mk_k_mat(double* K_ele, double* K, int n, int dof, int rows);
void mk_k_ele(char* argv[], double l, double* K, int rows);
void print_m(double* M, int rows);
int addr(int r, int c, int cols);

int main(int argc, char* argv[]) {
    cout << argc << endl;

    // Units are converted from those supplied to SI units
    double L = atof(argv[1]) * 0.001;
    const int n = atoi(argv[2]);
    const int N = n + 1;
    const int dof = 3;
    const int ele =  dof * 2;
    const int rows = N*dof;

    if(n%2 != 0){
      cout << "The number of elements must be even to ensure a node is present"
      << " below the concentrated force. Program quitting" << endl;
      return 0;
    }
    double l = L/n;

    cout << "Matrix size: " << rows << endl;

    double rho = 7850;
    double * K_ele = new double[ele*ele];
    double * K = new double[rows*rows];
    // mk_k_ele(argv, l, K_ele, ele);
    mk_k_ele(argv, l, K_ele, ele);
    print_m(K_ele, ele);
    cout << endl << endl;
    mk_k_mat(K_ele, K, n, dof, rows);
    print_m(K, rows);
}

void mk_k_mat(double* K_ele, double* K, int n, int dof, int rows){
    for (int a = 0; a < n; a++){
        for (int r = 1; r < dof*2+1; r++){
            for (int c = 1; c < dof*2+1; c++){
                K[addr(r + a * dof, c + a * dof, rows)] += K_ele[addr(r,c,dof*2)];
            }
        }
    }
}

void mk_k_ele(char* argv[], double l, double* K, int rows){
    // Units are converted from those supplied to SI units
    double A = atof(argv[3]) * pow(10, -6);
    double I = atof(argv[4]) * pow(10, -12);
    double E = atof(argv[5]) * pow(10, 6);
    // K will be stored in compressed form, so will be the elemental form
    K[addr(1,1,rows)] = A*E/l;
    K[addr(1,4,rows)] = -K[addr(1,1,rows)];
    K[addr(2,2,rows)] = 12 * E * I / pow(l, 3);
    K[addr(2,3,rows)] = 6 * E * I / pow(l, 2);
    K[addr(2,5,rows)] = -K[addr(2,2,rows)];
    K[addr(2,6,rows)] = K[addr(2,3,rows)];
    K[addr(3,3,rows)] = 4 * E * I / l;
    K[addr(3,5,rows)] = -K[addr(2,3,rows)];
    K[addr(3,6,rows)] = 2 * E * I / l;
    K[addr(4,4,rows)] = K[addr(1,1,rows)];
    K[addr(5,5,rows)] = K[addr(2,2,rows)];
    K[addr(5,6,rows)] = K[addr(3,5,rows)];
    K[addr(6,6,rows)] = K[addr(3,3,rows)];
    // Add symmetric terms
    for (int r = 2; r < rows + 1; r++){
        for (int c = 1; c < r; c++){
            K[addr(r,c,rows)] = K[addr(c,r,rows)];
        }
    }
}

void mk_k_ele_ONES(char* argv[], double l, double* K, int rows){
    // Make ones everywhere
    for (int r = 1; r < rows+1; r++){
        for (int c = 1; c < rows+1; c++){
            K[addr(r,c,rows)] = 1;
        }
    }
}

int addr(int r, int c, int cols){
    // Helper function to index from human-indices to c++ vector notation
    return (r-1)*cols+(c-1);
}



void print_m(double* M, int rows){
    for (int c = 0; c < rows; c++){
        for (int r = 0; r < rows; r++){
            cout << M[c*rows+r] << "\t";
        }
        cout << endl << endl;
    }
}
