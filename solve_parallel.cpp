#include <iostream>
#include <cstdlib>
#include <cmath>
#include <mpi.h>

using namespace std;

#include "helpers.h"
#include "mat_builder.h"
#include "lapack_c.h"


void solve_explicit_parallel(int argc, char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs, int bw){
    // Get from args the core variables and recalculate values
    double l = atof(argv[1]) * 0.001 / atoi(argv[2]);
    const double A = atof(argv[3]) * pow(10, -6);
    const int rows = 2 * bw + 1;
    const double T = atof(argv[7]);
    int iters = atof(argv[8]);
    double t_step = T / iters;
    const double rho = atof(argv[9]);

    int retval;
    retval = MPI_Init(&argc, &argv);
    cout << "Hello world" << endl;

    // can i see the variables outside?

    MPI_Finalize();


    // HERE THINGS CHANGE....

    // // F is the holder of the RHS terms as they are added and hold the u(t+1) solution at the end of each solve
    // // iteration, but will then have its pointer overwritten to point at the to-be-discareded u(t-1)
    // double* F = new double[eqs]();
    // double* u_past = new double[eqs]();
    // double* u_pres = new double[eqs]();
    //
    // // Memory allocation to default values needed by dgbsv
    // int* ipiv = new int[eqs]();
    // int info = 0;
    //
    // // The M_ref matrix is unmultiplied by its factor rho * A * l, so it is copied and multiplied at this stage.
    // // It gets also multiplied by 1 / ∆t^2 as all elements in the equation present this factor
    // const double t_fact =  rho * A * l / (t_step * t_step);
    // double* M = new double[eqs]();
    // F77NAME(dcopy) (eqs, M_ref, 1, M, 1);
    // F77NAME(dscal) (eqs, t_fact, M, 1);
    //
    // // Create matrix (K - 2 * M) and remove the top rows of zeros that are unnecesarily present for this solve routine
    // double* KM = new double[eqs * rows]();
    // mk_km_mat(KM, K_ref, M, eqs, bw);
    //
    // // Iterate x steps until T is reached
    // for (int iter = 1; iter < iters + 1; iter++){
    //     // t is necessary to calculate the iteration's force applied
    //     double t = iter * t_step;
    //
    //     // Overwrite the previous iteration's t-1 positions with the F vector
    //     F77NAME(dcopy) (eqs, F_ref, 1, F, 1);
    //     // and scale it by the iteration's time, but do not increase past t = 1.0s
    //     F77NAME(dscal) (eqs, min(t,1.0), F, 1);
    //
    //     // Multiply the t+0 position by the KM matrix and add to the force vector
    //     F77NAME(dgbmv)('t', eqs, eqs, bw, bw, -1, KM, rows, u_pres, 1, 1, F, 1);
    //
    //     // Multiply the t-1 position by M and add to the vector calculated above
    //     F77NAME(dgbmv)('N', eqs, eqs, 0, 0, -1, M, 1, u_past, 1, 1, F, 1);
    //
    //
    //     // F now holds the entire RHS and the equation can be solved as Ax = b
    //     int info = 0;
    //     F77NAME(dgbsv) (eqs, 0, 0, 1, M, 1, ipiv, F, eqs, &info);
    //     if (info){
    //         cout << "ERROR, dgbsv threw an exception: " << info << endl;
    //     }
    //
    //     // Instead of reallocating memory or inserting values in each vector, swap the pointers
    //     // to the arrays (which are of same size)
    //     shift_vec(u_past, u_pres, F);
    // }
    //
    // // As per exercise 1, write the solution to the F vector passed in by main()
    // F77NAME(dcopy) (eqs, u_pres, 1, F_ref, 1);
}
