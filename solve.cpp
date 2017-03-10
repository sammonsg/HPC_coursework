#include <iostream>
#include <cstdlib>

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
}

void solve_static(double* K, double* F, int eqs, int bw){
    double* y = new double[eqs];
    int info = 0;
    int* ipiv = new int[eqs]();
    int rows = 3 * bw + 1;
    F77NAME(dgbmv)('N', eqs, eqs, bw, bw, 1, K, rows, F, 1, 1, y, 1);
    F77NAME(dgbsv) (eqs, bw, bw, 1, K, rows, ipiv, F, eqs, &info);
    print_v(F, eqs);
}

void solve_explicit(char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs){
    double T = atof(argv[7]);
    int iters = atof(argv[8]);
    double t_step = 1.0 / iters;
    double rho = atof(argv[9]);

    for (int iter = 0; iter < iters; iter++){
        double t = iter * t_step;
        double * K = new double[eqs * 13]();
        double * F = new double[eqs]();
        double * M = new double[eqs]();
        // WORK IN PROGRESS

        // copy_n(K, eqs * 13, K_ref);
        // copy_n(F, eqs, F_ref);
        // copy_n(M, eqs, M_ref);
        // print_banded_m(K, 13, eqs);
        // print_banded_m(K_ref, 13, eqs);
        //////////////////
        cout << "The time is: " << t << endl;
    }
}
