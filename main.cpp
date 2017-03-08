#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "mat_builder.h"
#include "helpers.h"

int main(int argc, char* argv[]) {

    // Units are converted from those supplied to SI units
    double L = atof(argv[1]) * 0.001;
    const int n = atoi(argv[2]);
    const int N = n - 1;
    const int dof = 3;
    const int ele =  dof * 2;
    const int rows = N*dof;
    const int bw = 4;

    if(n%2 != 0){
      cout << "The number of elements must be even to ensure a node is present"
      << " below the concentrated force. Program quitting" << endl;
      return 0;
    }
    double l = L/n;

    cout << "Matrix size: " << rows << endl;

    double rho = 7850;
    double * K_ele = new double[ele*ele];
    double * K = new double[rows*(2*bw+1)];
    mk_banded_k_ele(argv, l, K_ele, ele, bw);
    // cout << endl << endl;
    // mk_banded_k_mat(K_ele, K, n, dof, rows, bw);
    // print_m(K, rows);
}

int full_m_solve(int argc, char* argv[]) {

    // Units are converted from those supplied to SI units
    double L = atof(argv[1]) * 0.001;
    const int n = atoi(argv[2]);
    const int N = n - 1;
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
    // print_m(K_ele, ele);
    cout << endl << endl;
    mk_k_mat(K_ele, K, n, dof, rows);
    print_m(K, rows);
    return 0;
}
