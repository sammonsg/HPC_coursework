#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "mat_builder.h"
#include "helpers.h"
#include "solve.h"

int main(int argc, char* argv[]) {
    // Get exercise number
    int ex = atoi(argv[6]);
    cout << "Exercise " << ex << endl;

    // Set output verbosity
    bool verbose = 0;

    // Units are converted from those supplied to SI units
    double L = atof(argv[1]) * 0.001;
    const int n = atoi(argv[2]);
    const int N = n - 1;
    const int dof = 3;
    const int ele =  dof * 2;
    const int cols = N*dof;
    const int bw = 4;
    const int rows = 3 * bw + 1;
    const int q_x = 0;
    const int q_y = 1000;
    const int F_centre = 1000;
    const double beta = 0.25;
    const double gamma = 0.5;

    if(n % 2 != 0){
      cout << "The number of elements must be even to ensure a node is present"
      << " below the concentrated force. Program quitting" << endl;
      return 0;
    }
    double l = L/n;

    double rho = 7850;

    // Allocate memory for K, F and M. There is negligible slowdown for the
    // exercises where M is not neeed.
    double * K = new double[cols*rows]();
    double * F = new double[cols]();
    // The M matrix is a diagonal matrix, so it can be expressed as a vector
    double * M = new double[cols]();
    mk_banded_k_mat(argv, l, K, N, dof, rows, bw);
    mk_F_mat(F, N, dof, q_x, q_y, F_centre, l);
    mk_M_mat(argv, M, N, l, dof);
    // mk_M_mat(char* argv[], double* M, int N, double l, int dof)


    if (verbose == true){
        cout << "Matrix rows: " << rows << "\tMatrix cols: " << cols << endl;
        print_banded_m(K, rows, dof*N);
        print_v(F, cols);
        print_v(M, cols);
    }
    if (ex == 1){
        solve_static(K, F, cols, bw);
        print_v(F, cols);

    }
    if (ex == 2){
        if (argc == 10){
            solve_explicit(argv, K, F, M, cols, bw);
        }
        else {
            cout << "Error, insufficient args supplied for exercise 2" << endl;
        }
    }

    if (ex == 2){
        if (argc == 10){
            solve_explicit(argv, K, F, M, cols, bw);
        }
        else {
            cout << "Error, insufficient args supplied for exercise 2" << endl;
        }
    }


    delete [] K;
    delete [] F;
    delete [] M;
    return 0;
}
