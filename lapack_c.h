// All extern calls are put into this file for reuse in other files


#define F77NAME(x) x##_
extern "C" {
    double F77NAME(dnrm2)(const int& N, const double *X, const int& incX);
    void F77NAME(dgbmv)(const char& trans, const int& m, const int& n, const int& kl, const int& ku,
                const double& alpha, const double* a, const int& lda, const double* x, const int& incx,
                const double& beta, double* y, const int& incy);
    void F77NAME(dgbsv)(const int& n, const int& kl, const int& ku, const int& nrhs, const double * A, const int& ldab,
                int * ipiv, double * B, const int& ldb, int* info);
    void F77NAME(dscal)(const int& n, const double& da, double * dx, const int& incx);
    void F77NAME(dcopy)(const int& n, double * dx, const int& incx, double * dy, const int& incy);
    void F77NAME(daxpy)(const int& n, const double& da, double * dx, const int& incx, double * dy, const int& incy);
    void F77NAME(dgemm)(const char& transa, const char& transb, const int& m, const int& n, const int& k,
                const double& alpha, const double* a, const int& lda, const double * b, const int& ldb,
                const double& beta, const double * c, const int& ldc);
    void F77NAME(pdgbsv)(const int& n, const int& kl, const int& ku,
                const int& nrhs, const double * A, const int& ja,
                const int* desca, int * ipiv, double * B, const int& ib,
                const int* descb, double* work, const int& lwork,
                int* info);
    void F77NAME(pdgbtrf) (const int& n, const int& bwl, const int& bwu,
            const double* a,     const int& ja,
            const int* desca,    int* ipiv,
            double* AF, const int& LAF, double* work, const int& LWORK,
            int& info);
    void F77NAME(pdgbtrs) (const char& trans, const int& n, const int& kl,
              const int &ku, const int& nrhs,   const double* a,
              const int& ja, const int* desca, const int* ipiv, double* b,
              const int& ib, const int* descb, double* AF, const int& LAF,
              double* work, const int& LWORK, int& info);
    void Cblacs_get(int, int, int*);
    void Cblacs_pinfo(int*, int*);
    void Cblacs_gridinit(int*, char*, int, int);
    void Cblacs_gridinfo(int, int*, int*, int*, int*);
    void Cblacs_exit(int);
    void Cblacs_gridexit(int);
}
