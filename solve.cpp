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
}

void solve_static(double* K, double* F, int eqs, int bw){
    double* y = new double[eqs];
    int info = 0;
    int* ipiv = new int[eqs]();
    int rows = 3 * bw + 1;
    F77NAME(dgbmv)('N', eqs, eqs, bw, bw, 1, K, rows, F, 1, 1, y, 1);
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
    const double M_fact = rho * A * l;

    // Iterate x steps until T is reached
    for (int iter = 0; iter < iters; iter++){
        double t = iter * t_step;

        // Clone
        // double * K = new double[eqs * rows]();
        // double * F = new double[eqs]();
        // double * M = new double[eqs]();
        //
        // copy_n(F_ref, eqs, F);
        // copy_n(M_ref, eqs, M);
        cout << "The time is: " << t << endl;
    }
}
