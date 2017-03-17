#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "helpers.h"
#include "mat_builder.h"
#include "lapack_c.h"

void u1_solve_routine(double* u0, double* v0, double* a0, double* F, double* M, double* Keff, double b_dt,
                      double b_dt2, double c_a0, int eqs, int K_rows, int* ipiv, int info);
void a1_solve_routine(double* u0, double* u1, double* v0, double* a0, double b_dt, double b_dt2, double c_a0, int eqs);
void v1_solve_routine(double* v0, double* a0, double* a1, double dt, double gmm, int eqs);

void solve_static(double* K, double* F, int eqs, int bw){
    int info = 0;
    int* ipiv = new int[eqs]();
    int rows = 3 * bw + 1;
    F77NAME(dgbsv) (eqs, bw, bw, 1, K, rows, ipiv, F, eqs, &info);
}

void solve_explicit(char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs, int bw){
    // Get from args the core variables and recalculate values
    double l = atof(argv[1]) * 0.001 / atoi(argv[2]);
    const double A = atof(argv[3]) * pow(10, -6);
    const int rows = 2 * bw + 1;
    const double T = atof(argv[7]);
    int iters = atof(argv[8]);
    double t_step = 1.0 / iters;
    const double rho = atof(argv[9]);

    // F is the holder of the RHS terms as they are added and hold the u(t+1) solution at the end of each solve
    // iteration, but will then have its pointer overwritten to point at the to-be-discareded u(t-1)
    double* F = new double[eqs]();
    double* u_past = new double[eqs]();
    double* u_pres = new double[eqs]();

    // Memory allocation to default values needed by dgbsv
    int* ipiv = new int[eqs]();
    int info = 0;

    // The M_ref matrix is unmultiplied by its factor rho * A * l, so it is copied and multiplied at this stage.
    // It gets also multiplied by 1 / ∆t^2 as all elements in the equation present this factor
    const double t_fact =  rho * A * l / (t_step * t_step);
    double* M = new double[eqs]();
    F77NAME(dcopy) (eqs, M_ref, 1, M, 1);
    F77NAME(dscal) (eqs, t_fact, M, 1);

    // Create matrix (K - 2 * M) and remove the top rows of zeros that are unnecesarily present for this solve routine
    double* KM = new double[eqs * rows]();
    mk_km_mat(KM, K_ref, M, eqs, bw);

    // Iterate x steps until T is reached
    for (int iter = 1; iter < iters + 1; iter++){
        // t is necessary to calculate the iteration's force applied
        double t = iter * t_step;

        // Overwrite the previous iteration's t-1 positions with the F vector
        F77NAME(dcopy) (eqs, F_ref, 1, F, 1);
        // and scale it by the iteration's force
        F77NAME(dscal) (eqs, t, F, 1);

        // Multiply the t+0 position by the KM matrix and add to the force vector
        F77NAME(dgbmv)('t', eqs, eqs, bw, bw, -1, KM, rows, u_pres, 1, 1, F, 1);

        // Multiply the t-1 position by M and add to the vector calculated above
        F77NAME(dgbmv)('N', eqs, eqs, 0, 0, -1, M, 1, u_past, 1, 1, F, 1);


        // F now holds the entire RHS and the equation can be solved as Ax = b
        int info = 0;
        F77NAME(dgbsv) (eqs, 0, 0, 1, M, 1, ipiv, F, eqs, &info);
        if (info){
            cout << "ERROR, dgbsv threw an exception: " << info << endl;
        }

        // Instead of reallocating memory or inserting values in each vector, swap the pointers
        // to the arrays (which are of same size)
        shift_vec(u_past, u_pres, F);
    }

    // As per exercise 1, write the solution to the F vector passed in by main()
    F77NAME(dcopy) (eqs, u_pres, 1, F_ref, 1);
}

void solve_implicit(char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs, int bw, double bt, double gmm){
    // Get from args the core variables and recalculate values
    double l = atof(argv[1]) * 0.001 / atoi(argv[2]);
    const double A = atof(argv[3]) * pow(10, -6);
    const int rows = 2 * bw + 1;
    const int K_rows = 3 * bw + 1;
    const double T = atof(argv[7]);
    int iters = atof(argv[8]);
    double dt = T / iters;
    const double rho = atof(argv[9]);
    const double M_fact =  rho * A * l;
    // Copy and scale M_ref by rho * A * l
    double* M = new double[eqs]();
    F77NAME(daxpy) (eqs, M_fact, M_ref, 1, M, 1);

    print_banded_m(K_ref, K_rows, eqs);

    // Calculate common coefficients
    double b_dt = 1 / (bt * dt);
    double b_dt2 = b_dt / dt;
    double c_a0 = (1 / (2 * bt)) - 1;



    // Create Keff matrix and remove the top rows of zeros that are unnecesarily present for this solve routine
    double* Keff = new double[eqs * K_rows]();
    mk_keff_mat(Keff, K_ref, M, b_dt2, eqs, bw);

    // Position vectors, u0 = u(t), u1 = u(t+1)
    double* u0 = new double[eqs]();
    double* u1 = new double[eqs]();

    // Velocity vectors, v0 = v(t), v1 = v(t+1)
    double* v0 = new double[eqs]();
    double* v1 = new double[eqs]();

    // Acceleration vectors, a0 = a(t), a1 = a(t+1)
    double* a0 = new double[eqs]();
    double* a1 = new double[eqs]();

    double* temp = new double[eqs]();
    int* ipiv = new int[eqs]();
    int info = 0;

    for (int iter = 0; iter < iters + 1; iter++){
        double t = iter * dt;
        F77NAME(dcopy) (eqs, F_ref, 1, temp, 1);
        // Scale by the time. Use min(t, 1) to get the effect of linearly increasing load until t=1 and then hold force
        F77NAME(dscal) (eqs, min(t+dt,1.0), temp, 1);

        // Calculate u(t+1) into temp and copy into u1
        u1_solve_routine(u0, v0, a0, temp, M, Keff, b_dt, b_dt2, c_a0, eqs, K_rows, ipiv, info);
        F77NAME(dcopy) (eqs, temp, 1, u1, 1);

        // Calculate a(t+1) into temp and copy into a1
        a1_solve_routine(u0, temp, v0, a0, b_dt, b_dt2, c_a0, eqs);
        F77NAME(dcopy) (eqs, temp, 1, a1, 1);

        // Calculate v(t+1) into temp and copy into v1
        v1_solve_routine(v0, a0, temp, dt, gmm, eqs);
        F77NAME(dcopy) (eqs, temp, 1, v1, 1);


        // Swap pointers to get f(t+1) into f(t) position
        shift_vec(u0, u1);
        shift_vec(v0, v1);
        shift_vec(a0, a1);

        if (info){
            break;
        }
    }
    F77NAME(dcopy) (eqs, u0, 1, F_ref, 1);
    print_pos_v(F_ref, eqs);
}

void u1_solve_routine(double* u0, double* v0, double* a0, double* F, double* M, double* Keff, double b_dt,
                      double b_dt2, double c_a0, int eqs, int K_rows, int* ipiv, int info){
    // Create new vector, copy a0 scaled by (1/2β)-1
    double* tmp = new double[eqs * K_rows]();
    F77NAME(daxpy) (eqs, c_a0, a0, 1, tmp, 1);

    // Add the v0 term to tmp, v0 scaled by 1/(β * dt)
    F77NAME(daxpy) (eqs, b_dt, v0, 1, tmp, 1);

    // Add the u0 term to tmp, u0 scaled by 1/(β * dt * dt)
    F77NAME(daxpy) (eqs, b_dt2, u0, 1, tmp, 1);

    // print_v(F, eqs);

    // Get F as result to Ax + y, then dump tmp
    F77NAME(dgbmv) ('n', eqs, eqs, 0, 0, 1, M, 1, tmp, 1, 1, F, 1);


    // print_banded_m(Keff, K_rows, eqs);
    // print_v(F, eqs);
    F77NAME(dcopy) (eqs * K_rows, Keff, 1, tmp, 1);
    // Solve the system into F
    F77NAME(dgbsv) (eqs, 4, 4, 1, tmp, K_rows, ipiv, F, eqs, &info);
    delete [] tmp;

    // print_v(F, eqs);
    // print_pos_v(F, eqs);

    if (info){
        cout << "An error occurred in dgbsv during solve of u1" << endl;
    }
}
void a1_solve_routine(double* u0, double* u1, double* v0, double* a0, double b_dt, double b_dt2, double c_a0, int eqs){
    // Subtract u0 from u1, and scale by its factor
    // Note: u1 within this function is a temporary variable, there is no overwriting of u1 outside this function
    F77NAME(daxpy) (eqs, -1, u0, 1, u1, 1);
    F77NAME(dscal) (eqs, b_dt2, u1, 1);

    // Add (u1 - u0) to the velocity term
    F77NAME(daxpy) (eqs, -b_dt, v0, 1, u1, 1);

    // Add the accelertion term
    F77NAME(daxpy) (eqs, -c_a0, a0, 1, u1, 1);
}
void v1_solve_routine(double* v0, double* a0, double* a1, double dt, double gmm, int eqs){
    // Solve routine onto a1, which is a temporary variable outside of this scope
    F77NAME(dscal) (eqs, dt * gmm, a1, 1);

    // Add the a0 term
    F77NAME(daxpy) (eqs, dt*(1-gmm), a0, 1, a1, 1);

    // Add the v0 term
    F77NAME(daxpy) (eqs, dt*(1-gmm), v0, 1, a1, 1);
}
