#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "helpers.h"
#include "mat_builder.h"
#include "lapack_c.h"
#include "solve.h"
#include <mpi.h>

void u1_solve_parallel(double* u0, double* v0, double* a0, double* F, double* M, double* Keff, double b_dt,
                       double b_dt2, double c_a0, int eqs, int K_rows, int* ipiv, int info);
void fill_descriptors(int n, int nb, int bw, int ctx, int* desca, int* descb);

void get_solve_domain(int eqs, int rank, int cores, int* begin, int* end){
    int std_size = eqs / cores + 1;
    int extra_cols = eqs % cores;
    *begin = rank * std_size - max((rank - extra_cols),0);
    *end = *begin + std_size - 1;
    if (rank >= extra_cols){
        *end =  *end -1;
    }
}

void exchange_boundaries(double* F, int solvespace, bool RHS_edge, bool LHS_edge, int rank, int iter){
    if(!RHS_edge){
        MPI_Send(F+(solvespace-8), 4, MPI_DOUBLE, rank + 1, iter, MPI_COMM_WORLD);
        MPI_Recv(F+(solvespace-4), 4, MPI_DOUBLE, rank + 1, iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(!LHS_edge){
        MPI_Send(F+4, 4, MPI_DOUBLE, rank - 1, iter, MPI_COMM_WORLD);
        MPI_Recv(F, 4, MPI_DOUBLE, rank - 1, iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void gather_solution(double* F_ref, int eqs, double* u_pres, int rank, int cores, int solvespace){
    if(rank == 0){
        F77NAME(dcopy) (solvespace-4, u_pres, 1, F_ref, 1);
        for (int core = 1; core < cores; core++){
            int begin = 0;
            int end = 0;
            get_solve_domain(eqs, core, cores, &begin, &end);
            int len = end - begin + 1;
            MPI_Recv(F_ref + begin, len, MPI_DOUBLE, core, core, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //Receive data from all cores
    }
    if(rank == cores - 1){
        MPI_Send(u_pres+4, solvespace - 4, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
    }
    else{
        MPI_Send(u_pres+4, solvespace - 8, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
    }
}

void solve_explicit_parallel(char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs, int bw,
                            int rank, int cores){

    // Get from args the core variables and recalculate values
    double l = atof(argv[1]) * 0.001 / atoi(argv[2]);
    const double A = atof(argv[3]) * pow(10, -6);
    const int rows = 2 * bw + 1;
    const int K_rows = rows + bw;
    const double T = atof(argv[7]);
    int iters = atof(argv[8]);
    double t_step = T / iters;
    const double rho = atof(argv[9]);

    // Calculate the solve-domain of each rank
    int begin = 0;
    int end = 0;
    get_solve_domain(eqs, rank, cores, &begin, &end);


    //Get the breadth of the domain of influence and make flags for edges
    bool RHS_edge = 0;
    bool LHS_edge = 0;
    if (rank == 0) {LHS_edge = 1;}
    if (rank == cores-1) {RHS_edge = 1;}
    int lhs_solvespace = max(begin - 4 + 4 * LHS_edge, 0);
    int rhs_solvespace = min(end + 4 - 4 * RHS_edge, eqs-1);
    int solvespace = rhs_solvespace - lhs_solvespace + 1;

    // cout << rank << " has been allocated " << begin << " to " << end << endl;
    // cout << rank << " has influence across " << lhs_solvespace << " to " << rhs_solvespace << endl;

    // Allocate memory for rank's solve space
    double* rank_K = new double[solvespace * K_rows]();
    double* rank_F = new double[solvespace]();
    double* rank_M = new double[solvespace]();
    mk_truncated_mat(rank_K, K_ref, K_rows, lhs_solvespace, rhs_solvespace, 0);
    mk_truncated_v(rank_F, F_ref, lhs_solvespace, rhs_solvespace);
    mk_truncated_v(rank_M, M_ref, lhs_solvespace, rhs_solvespace);



    // F is the holder of the RHS terms as they are added and hold the u(t+1) solution at the end of each solve
    // iteration, but will then have its pointer overwritten to point at the to-be-discareded u(t-1)
    double* F = new double[solvespace]();
    double* u_past = new double[solvespace]();
    double* u_pres = new double[solvespace]();

    // The M_ref matrix is unmultiplied by its factor rho * A * l, so it is copied and multiplied at this stage.
    // It gets also multiplied by 1 / ∆t^2 as all elements in the equation present this factor
    const double t_fact =  rho * A * l / (t_step * t_step);
    double* M = new double[solvespace]();
    double* M_inv = new double[solvespace]();
    F77NAME(dcopy) (solvespace, rank_M, 1, M, 1);
    F77NAME(dscal) (solvespace, t_fact, M, 1);
    invert_v(M_inv, M, solvespace);

    // if(rank == 2){
    //     print_v(rank_F, solvespace);
    // }

    // Create matrix (K - 2 * M) and remove the top rows of zeros that are unnecesarily present for this solve routine
    double* KM = new double[solvespace * rows]();
    mk_km_mat(KM, rank_K, M, solvespace, bw);

    for (int iter = 1; iter < iters + 1; iter++){
        // Get the simulation time, to scale the force accordingly
        double t = iter * t_step;

        // Run encapsulated solver. It was encapsulated for use by the parallel explicit solver.
        solve_iter(u_past, u_pres, F, rank_F, KM, M, M_inv, t, rows, solvespace, bw);

        // Broadcast solution to neighbors and receive updated values
        exchange_boundaries(F, solvespace, RHS_edge, LHS_edge, rank, iter);

        // Swap out pointer locations
        shift_vec(u_past, u_pres, F);
    }

    // Gather the solution on core 0
    gather_solution(F_ref, eqs, u_pres, rank, cores, solvespace);

    // print_v(u_pres+4, 4);
    delete [] rank_K;
    delete [] rank_F;
    delete [] rank_M;
}

void solve_implicit_parallel(char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs, int bw,
                            double bt, double gmm, int rank, int cores){

    // Get from args the core variables and recalculate values
    double l = atof(argv[1]) * 0.001 / atoi(argv[2]);
    const double A = atof(argv[3]) * pow(10, -6);
    const int rows = 4 * bw + 1;
    const int K_rows = rows - bw;
    const double T = atof(argv[7]);
    int iters = atof(argv[8]);
    double dt = T / iters;
    const double rho = atof(argv[9]);
    int solvespace = (eqs) / cores;

    // Adjust size if domain is not divisible equally amongst cores
    if (eqs % cores) { solvespace++; }
    const int lhs_pt = rank * solvespace;
    const int rhs_pt = (rank + 1) * solvespace - 1;
    cout << "Core " << rank << " has been allocated " << lhs_pt << " to ";
    cout << rhs_pt << " for a total of " << solvespace << " in domain of size " << eqs << endl;
    int extra_cols = solvespace * cores - eqs;

    // Calculate common coefficients
    double b_dt = 1 / (bt * dt);
    double b_dt2 = b_dt / dt;
    double c_a0 = (1 / (2 * bt)) - 1;

    // Make local elements from full matrices
    double* rank_K = new double[solvespace * rows]();
    double* rank_F = new double[solvespace]();
    double* rank_M = new double[solvespace]();
    double* Keff = new double[solvespace * rows]();

    // Use of minimum prevents from copying beyond size of reference matrix
    mk_truncated_mat(rank_K, K_ref, rows, lhs_pt, min(rhs_pt,eqs-1), -4);
    mk_truncated_v(rank_F, F_ref, lhs_pt, min(rhs_pt,eqs-1));
    mk_truncated_v(rank_M, M_ref, lhs_pt, min(rhs_pt,eqs-1));

    // Scale rank_M by rho * A * l
    F77NAME(dscal) (solvespace, rho * A * l, rank_M, 1);
    mk_keff_mat(Keff, rank_K, rank_M, b_dt2, solvespace, bw, 4);

    // Correct for extended domain on last core
    if(rank == cores - 1){
        for (int c = solvespace - extra_cols; c < solvespace; c++){
            Keff[rows - bw - 1 + c * rows] = 1;
            rank_M[c] = 1;
        }
    }

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


        F77NAME(dcopy) (solvespace, F_ref, 1, temp, 1);
        // Scale by the time. Use min(t, 1) to get the effect of linearly increasing load until t=1 and then hold force
        F77NAME(dscal) (solvespace, min(t+dt,1.00), temp, 1);

        // Calculate u(t+1) into temp and copy into u1
        u1_solve_parallel(u0, v0, a0, temp, rank_M, Keff, b_dt, b_dt2, c_a0, solvespace, K_rows, ipiv, info);

        F77NAME(dcopy) (solvespace, temp, 1, u1, 1);

        // Calculate a(t+1) into temp and copy into a1
        a1_solve_routine(u0, temp, v0, a0, b_dt, b_dt2, c_a0, solvespace);
        F77NAME(dcopy) (solvespace, temp, 1, a1, 1);

        // Calculate v(t+1) into temp and copy into v1
        v1_solve_routine(v0, a0, temp, dt, gmm, solvespace);
        F77NAME(dcopy) (solvespace, temp, 1, v1, 1);


        // Swap pointers to get f(t+1) into f(t) position
        shift_vec(u0, u1);
        shift_vec(v0, v1);
        shift_vec(a0, a1);

        if (info){
            break;
        }
    }
    // F77NAME(dcopy) (solvespace, u0, 1, F_ref, 1);

}

void u1_solve_parallel(double* u0, double* v0, double* a0, double* F, double* M, double* Keff, double b_dt,
                       double b_dt2, double c_a0, int eqs, int K_rows, int* ipiv, int info){

    // Create new vector, copy a0 scaled by (1/2β)-1
    double* tmp = new double[eqs * K_rows]();
    F77NAME(daxpy) (eqs, c_a0, a0, 1, tmp, 1);

    // Add the v0 term to tmp, v0 scaled by 1/(β * dt)
    F77NAME(daxpy) (eqs, b_dt, v0, 1, tmp, 1);

    // Add the u0 term to tmp, u0 scaled by 1/(β * dt * dt)
    F77NAME(daxpy) (eqs, b_dt2, u0, 1, tmp, 1);

    // Get F as result to Ax + y, then dump tmp
    F77NAME(dgbmv) ('n', eqs, eqs, 0, 0, 1, M, 1, tmp, 1, 1, F, 1);

    F77NAME(dcopy) (eqs * K_rows, Keff, 1, tmp, 1);
    // Solve the system into F
    // F77NAME(pdgbsv) (eqs, 4, 4, 1, tmp, 1, desca, ipiv, F, eqs, &info);
    delete [] tmp;

    if (info){
        cout << "An error occurred in pdgbsv during PDGBSV: " << info << endl;
    }
}

void fill_descriptors(int n, int nb, int rows, int ctx, int* desca, int* descb){
    // Create descriptors for matrix and RHS vector storage
    desca[0] = 501;         // Type is a banded matrix 1-by-P
    desca[1] = ctx;         // Context
    desca[2] = n;           // Problem size
    desca[3] = nb;          // Blocking
    desca[4] = 0;           // Process row/column
    desca[5] = rows;         // Local size
    desca[6] = 0;           // Reserved

    descb[0] = 502;         // Type is a banded matrix P-by-1 (RHS)
    descb[1] = ctx;         // Context
    descb[2] = n;           // Problem size
    descb[3] = nb;          // Blocking
    descb[4] = 0;           // Process row/column
    descb[5] = nb;          // Local size
    descb[6] = 0;           // Reserved
}
