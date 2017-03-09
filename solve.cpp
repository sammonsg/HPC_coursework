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
