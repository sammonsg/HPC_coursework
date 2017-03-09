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
    // std::cout << " I LIVE " << std::endl;
    int info = 0;
    int* ipiv = new int[eqs]();
    int rows = 3 * bw + 1;
    F77NAME(dgbmv)('N', eqs, eqs, bw, bw, 1, K, rows, F, 1, 1, y, 1);
    F77NAME(dgbsv) (eqs, bw, bw, 1, K, rows, ipiv, F, eqs, &info);
    std::cout << " I have completed " << std::endl;
    print_v(F, eqs);
}
