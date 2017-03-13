#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "helpers.h"

#define F77NAME(x) x##_
extern "C" {
    double F77NAME(dnrm2)(const int& N, const double *X, const int& incX);
    void F77NAME(daxpy)(const int& N, const double& alpha, const double *X,
                         const int& incX, double *Y, const int& incY);
    void F77NAME(dgbmv)(const char& trans, const int& m, const int& n,
                const int& kl, const int& ku,
                const double& alpha, const double* a, const int& lda,
                const double* x, const int& incx, const double& beta,
                double* y, const int& incy);
    void F77NAME(dgbsv)(const int& n, const int& kl, const int& ku,
                const int& nrhs, const double * A, const int& ldab,
                int * ipiv, double * B, const int& ldb, int* info);
    void F77NAME(dscal)(const int& n, const double& da, double * dx,
                const int& incx);
    void F77NAME(dcopy)(const int& n, double * dx, const int& incx,
                double * dy, const int& incy);
}

void solve_static(double* K, double* F, int eqs, int bw){
    int info = 0;
    int* ipiv = new int[eqs]();
    int rows = 3 * bw + 1;
    F77NAME(dgbsv) (eqs, bw, bw, 1, K, rows, ipiv, F, eqs, &info);
}

void solve_explicit(char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs, int bw){
    // Recover from args the core variables
    double l = atof(argv[1]) * 0.001 / atoi(argv[2]);
    const double A = atof(argv[3]) * pow(10, -6);
    const int rows = 3 * bw + 1;
    const double T = atof(argv[7]);
    int iters = atof(argv[8]);
    double t_step = 1.0 / iters;
    const double rho = atof(argv[9]);
    // cout << " RHO A L = " << M_fact << endl;
    const double t_fact =  rho * A * l / (t_step * t_step);
    cout << "Time_factor = " << t_fact << endl;
    double* u_past = new double[eqs]();
    double* u_pres = new double[eqs]();
    double* u_futu = new double[eqs]();
    double* M = new double[eqs]();
    F77NAME(dcopy) (eqs, M_ref, 1, M, 1);
    F77NAME(dscal) (eqs, t_fact, M, 1);

    // double* X = new double[eqs * rows]();
    // double* M = new double[eqs * rows]();

    // clone_vector(M_ref, M, eqs * rows);
    // print_v(M_ref, eqs);
    print_v(M, eqs);

    // m_diag_add(M, M_ref, t_fact, 1, 1);
    print_v(M_ref, eqs);

    // clone_vector(K_ref, X, eqs * rows);


    // print_v(M_ref, eqs);
    // print_banded_m(X, rows, eqs);

    // m_diag_add(X, M_ref, -2 * t_fact, eqs, bw);
    // print_banded_m(X, rows, eqs);
    // print_banded_m(K_ref, rows, eqs);





    // // Iterate x steps until T is reached
    // for (int iter = 0; iter < 4; iter++){
    //     clone_vector(F_ref, u_futu, eqs);
    //     // print_v(u_futu, eqs);
    //
    //     F77NAME(dgbmv)('N', eqs, eqs, bw, bw, -1, X, rows, u_pres, 1, 1, u_pres, 1);
    //     F77NAME(dgbmv)('N', eqs, eqs, bw, bw, -1, X, rows, u_past, 1, 1, u_futu, 1);
    //
    //     int info = 0;
    //     int* ipiv = new int[eqs]();
    //     print_v(u_futu, eqs);
    //
    //     F77NAME(dgbsv) (eqs, bw, bw, 1, M_ref, rows, ipiv, u_futu, eqs, &info);
    //
    //     if (info){
    //         cout << "ERROR, dgbsv threw an exception: " << info << endl;
    //     }
    //
    //         // print_v(u_futu, eqs);
    //
    //     shift_vec(u_past, u_pres, u_futu);
    //
    //
    //
    //     double t = iter * t_step;
    //
    //
    //     // Clone
    //
    //     // double * K = new double[eqs * rows]();
    //     // double * F = new double[eqs]();
    //     // double * M = new double[eqs]();
    //     //
    //     // copy_n(F_ref, eqs, F);
    //     // copy_n(M_ref, eqs, M);
    //
    //     // cout << "The time is: " << t << endl;
    // }
}
