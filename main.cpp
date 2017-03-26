#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "mat_builder.h"
#include "helpers.h"
#include "solve.h"
#include "solve_parallel.h"
#include <mpi.h>
#include "io.h"


int main(int argc, char* argv[]) {


    int retval;
    int rank = 0;
    int cores = 1;
    int ex = atoi(argv[6]);

    if (ex > 3){
        retval = MPI_Init(&argc, &argv);
        int retval;
        MPI_Comm_rank(MPI_COMM_WORLD,&rank ) ; // Get core number
        MPI_Comm_size(MPI_COMM_WORLD,&cores ) ; // Get core count
    }
    if (rank == 0){
    // Get exercise number
        cout << "Exercise " << ex << endl;
        printHeader();
    }

    // Set output verbosity
    bool show_matrices = 0;

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

    printInfo(n, cols, rank, cores);

    if(n % 2 != 0){
      cout << "The number of elements must be even to ensure a node is present"
      << " below the concentrated force. Program quitting" << endl;
      return 0;
    }
    double l = L/n;

    double rho = 7850;

    // Allocate memory for K, F and M. There is negligible slowdown for the exercises where M is not neeed.
    double * K = new double[cols*rows]();
    double * F = new double[cols]();
    // The M matrix is a diagonal matrix, so it can be expressed as a vector
    double * M = new double[cols]();
    mk_banded_k_mat(argv, l, K, N, dof, rows, bw);
    mk_F_mat(F, N, dof, q_x, q_y, F_centre, l);
    mk_M_mat(argv, M, N, l, dof);

    // Flag for showing solution on terminal
    bool print_solution = true;

    if (show_matrices == true){
        cout << "Matrix rows: " << rows << "\tMatrix cols: " << cols << endl;
        print_banded_m(K, rows, dof*N);
        print_v(F, cols);
        print_v(M, cols);
    }
    if (ex == 1){
        solve_static(K, F, cols, bw);
        if (print_solution) {
            print_pos_v(F, cols);
        }
    }
    if (ex == 2){
        if (argc == 10){

            // Turn this flag to True to get the load_rate study
            bool load_rate_study = false;

            if (load_rate_study) {
                int iters = atoi(argv[8]);
                int runs = iters/50 + 1;
                double* minmax = new double[iters]();
                for (int i = 0; i < runs - 1; i++){
                    double T_load = (double)i*50 / ((double)iters / 5.0);
                    // cout << i << "\t" << T_load << endl;
                    minmax[i*3] = i;
                    minmax[i*3+1] = T_load;
                    minmax[i*3+2] = solve_explicit(argv, K, F, M, cols, bw, T_load);
                    cout << "Iter: " << i << "\tT_load: " << T_load << "\tOscillations: " << minmax[3*i+2] << endl;
                }
                write_v("minmax_exp", minmax, runs, 3, false);
            }
            solve_explicit(argv, K, F, M, cols, bw);
            if (print_solution) {
                print_pos_v(F, cols);
            }
        }
        else {
            cout << "Error, wrong number of args supplied for exercise 2" << endl;
        }
    }

    if (ex == 3){

        if (argc == 10){
            bool load_rate_study = false;

            // Logic block for load rates
            if (load_rate_study) {
                int iters = atoi(argv[8]);
                int runs = iters/10 + 1;
                double* minmax = new double[iters]();
                for (int i = 0; i < runs - 1; i++){
                    double T_load = (double)i*10 / ((double)iters / 5.0);
                    // cout << i << "\t" << T_load << endl;
                    minmax[i*3] = i;
                    minmax[i*3+1] = T_load;
                    minmax[i*3+2] = solve_implicit(argv, K, F, M, cols, bw, beta, gamma, T_load);
                    cout << "Iter: " << i << "\tT_load: " << T_load << "\tOscillations: " << minmax[3*i+2] << endl;
                }
                write_v("minmax_imp", minmax, runs, 3, false);
            }

            // In any case, solve the standard solve, for verification purposes
            solve_implicit(argv, K, F, M, cols, bw, beta, gamma);
            if (print_solution){
                print_pos_v(F, cols);
            }
        }
        else {
            cout << "Error, wrong number of args supplied for exercise 3" << endl;
        }
    }
    if (ex == 4){
        if (argc != 10){
            cout << "Error, wrong number of args supplied for exercise 4" << endl;
            exit (EXIT_FAILURE);
        }
        solve_explicit_parallel(argv, K, F, M, cols, bw, rank, cores);
        if (rank == 0){
            if (print_solution){
                print_pos_v(F, cols);
            }
        }
    }
    if (ex == 5){
        if (argc != 10){
            cout << "Error, wrong number of args supplied for exercise 4" << endl;
            exit (EXIT_FAILURE);
        }
        solve_implicit_parallel(argv, K, F, M, cols, bw, beta, gamma, rank, cores);
        if (rank == 0){
            if (print_solution){
                print_pos_v(F, cols);
            }
        }
    }



    delete [] K;
    delete [] F;
    delete [] M;

    if (ex > 3){
        MPI_Finalize();
    }


    return 0;
}
