#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "helpers.h"
#include "mat_builder.h"
#include "lapack_c.h"
#include "solve.h"
#include <mpi.h>

void u1_solve_parallel(double* u0, double* v0, double* a0, double* F, double* M, double* Keff, double* tmp, double* wk,
                       double lwork, double b_dt, double b_dt2, double c_a0, int local_eqs, int eqs, int K_rows,
                       int bw, int* desca, int* descb, int rank, int* ipiv, int info);
void initialize_cblacs(int n, int nb, int bw, int ctx, int* desca, int* descb);

void get_solve_domain(int eqs, int rank, int cores, int* begin, int* end){
    int std_size = eqs / cores + 1;
    int extra_cols = eqs % cores;
    *begin = rank * std_size - max((rank - extra_cols),0);
    *end = *begin + std_size - 1;
    if (rank >= extra_cols){
        *end =  *end -1;
    }
}

void exchange_boundaries(double* F, int local_eqs, int bw, bool RHS_edge, bool LHS_edge, int rank, int iter){
    if(!RHS_edge){
        MPI_Send(F+(local_eqs-2*bw), bw, MPI_DOUBLE, rank + 1, iter, MPI_COMM_WORLD);
        MPI_Recv(F+(local_eqs-bw), bw, MPI_DOUBLE, rank + 1, iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(!LHS_edge){
        MPI_Send(F+bw, bw, MPI_DOUBLE, rank - 1, iter, MPI_COMM_WORLD);
        MPI_Recv(F, bw, MPI_DOUBLE, rank - 1, iter, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

void gather_solution(double* F_ref, int eqs, double* u_pres, int rank, int cores, int local_eqs, int bw, int blankCols = 0){
    if(rank == 0){
        F77NAME(dcopy) (local_eqs-bw, u_pres, 1, F_ref, 1);
        for (int core = 1; core < cores; core++){
            int begin;
            int end;
            if (blankCols == 0){
                get_solve_domain(eqs, core, cores, &begin, &end);
            }
            else{
                int true_solvespace = eqs / cores;
                if (eqs % cores) { true_solvespace++; }
                begin = core * true_solvespace;
                end = (core + 1) * true_solvespace - 1;
            }

            // Correction for last element
            int len = end - begin + 1;
            if (core == cores - 1 && blankCols){
                len -= blankCols;
            }
            // Correction for last element for
            MPI_Recv(F_ref + begin, len, MPI_DOUBLE, core, core, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        //Receive data from all cores
    }
    else if(rank == cores - 1){
        MPI_Send(u_pres+bw, local_eqs - bw - blankCols, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
    }
    else{
        MPI_Send(u_pres+bw, local_eqs - 2*bw, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

// gather_implicit_solution(F_ref, eqs, u0, rank, cores, local_eqs){
//}

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
    int local_eqs = rhs_solvespace - lhs_solvespace + 1;

    // cout << rank << " has been allocated " << begin << " to " << end << endl;
    // cout << rank << " has influence across " << lhs_solvespace << " to " << rhs_solvespace << endl;

    // Allocate memory for rank's solve space
    double* rank_K = new double[local_eqs * K_rows]();
    double* rank_F = new double[local_eqs]();
    double* rank_M = new double[local_eqs]();
    mk_truncated_mat(rank_K, K_ref, K_rows, lhs_solvespace, rhs_solvespace, 0);
    mk_truncated_v(rank_F, F_ref, lhs_solvespace, rhs_solvespace);
    mk_truncated_v(rank_M, M_ref, lhs_solvespace, rhs_solvespace);



    // F is the holder of the RHS terms as they are added and hold the u(t+1) solution at the end of each solve
    // iteration, but will then have its pointer overwritten to point at the to-be-discareded u(t-1)
    double* F = new double[local_eqs]();
    double* u_past = new double[local_eqs]();
    double* u_pres = new double[local_eqs]();

    // The M_ref matrix is unmultiplied by its factor rho * A * l, so it is copied and multiplied at this stage.
    // It gets also multiplied by 1 / ∆t^2 as all elements in the equation present this factor
    const double t_fact =  rho * A * l / (t_step * t_step);
    double* M = new double[local_eqs]();
    double* M_inv = new double[local_eqs]();
    F77NAME(dcopy) (local_eqs, rank_M, 1, M, 1);
    F77NAME(dscal) (local_eqs, t_fact, M, 1);
    invert_v(M_inv, M, local_eqs);



    // Create matrix (K - 2 * M) and remove the top rows of zeros that are unnecesarily present for this solve routine
    double* KM = new double[local_eqs * rows]();
    mk_km_mat(KM, rank_K, M, local_eqs, bw);

    for (int iter = 1; iter < iters + 1; iter++){
        // Get the simulation time, to scale the force accordingly
        double t = iter * t_step;

        // Run encapsulated solver. It was encapsulated for use by the parallel explicit solver.
        solve_iter(u_past, u_pres, F, rank_F, KM, M, M_inv, t, 1, rows, local_eqs, bw);

        // Broadcast solution to neighbors and receive updated values
        exchange_boundaries(F, local_eqs, bw, RHS_edge, LHS_edge, rank, iter);

        // Swap out pointer locations
        shift_vec(u_past, u_pres, F);


    }

    // Gather the solution on core 0
    gather_solution(F_ref, eqs, u_pres, rank, cores, local_eqs, bw);

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
    int local_eqs = eqs / cores;

    // Adjust size if domain is not divisible equally amongst cores
    if (eqs % cores) { local_eqs++; }
    const int lhs_pt = rank * local_eqs;
    const int rhs_pt = (rank + 1) * local_eqs - 1;
    int extra_cols = local_eqs * cores - eqs;

    // Calculate common coefficients
    double b_dt = 1 / (bt * dt);
    double b_dt2 = b_dt / dt;
    double c_a0 = (1 / (2 * bt)) - 1;

    // Make local elements from full matrices
    double* rank_K = new double[local_eqs * rows]();
    double* rank_F = new double[local_eqs]();
    double* rank_M = new double[local_eqs]();
    double* Keff = new double[local_eqs * rows]();

    // Use of minimum prevents from copying beyond size of reference matrix
    mk_truncated_mat(rank_K, K_ref, rows, lhs_pt, min(rhs_pt,eqs-1), -4);
    mk_truncated_v(rank_F, F_ref, lhs_pt, min(rhs_pt,eqs-1));
    mk_truncated_v(rank_M, M_ref, lhs_pt, min(rhs_pt,eqs-1));

    // Scale rank_M by rho * A * l
    F77NAME(dscal) (local_eqs, rho * A * l, rank_M, 1);
    mk_keff_mat(Keff, rank_K, rank_M, b_dt2, local_eqs, bw, 4);

    // Correct for extended domain on last core
    if(rank == cores - 1){
        for (int c = local_eqs - extra_cols; c < local_eqs; c++){
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

    // Carrier variable for use in solver
    double* temp = new double[eqs]();
    double* temp2 = new double[eqs*rows]();

    // Initialize Cblacs

    // Initialize scaLAPACK elements
    // Code snippet from pbsv.cpp as per Lecture 10 of HPC course at Imperial College of London
    int nrow = 1;
    int ncol = cores;
    char order = 'C';
    int ctx;
    int mype;
    int npe;
    int myrow;
    int mycol;
    Cblacs_pinfo(&mype, &npe);
    Cblacs_get( 0, 0, &ctx );
    Cblacs_gridinit( &ctx, &order, 1, npe );
    Cblacs_gridinfo( ctx, &nrow, &ncol, &myrow, &mycol);

    // Create descriptors for matrix and RHS vector storage
    int desca[7];
    desca[0] = 501;         // Type is a banded matrix 1-by-P
    desca[1] = ctx;         // Context
    desca[2] = eqs;           // Problem size
    desca[3] = local_eqs;          // Blocking
    desca[4] = 0;           // Process row/column
    desca[5] = rows;         // Local size
    desca[6] = 0;           // Reserved

    int descb[7];
    descb[0] = 502;         // Type is a banded matrix P-by-1 (RHS)
    descb[1] = ctx;         // Context
    descb[2] = eqs;           // Problem size
    descb[3] = local_eqs;          // Blocking
    descb[4] = 0;           // Process row/column
    descb[5] = local_eqs;          // Local size
    descb[6] = 0;           // Reserved

    const int lwork = 2*(local_eqs+bw)*(2*bw)+(12*bw)*(3*bw) + max(1*(local_eqs+6*bw), 1);
    int* ipiv = new int[eqs]();
    double* wk = new double[lwork];
    int info = 0;

    for (int iter = 1; iter < iters + 1; ++iter){
        double t = iter * dt;
        // cout << "t = " << t << endl;
        F77NAME(dcopy) (local_eqs, rank_F, 1, temp, 1);
        // Scale by the time. Use min(t, 1) to get the effect of linearly increasing load until t=1 and then hold force
        F77NAME(dscal) (local_eqs, min(t+dt,1.00), temp, 1);

        // Calculate u(t+1) into temp and copy into u1
        u1_solve_parallel(u0, v0, a0, temp, rank_M, Keff, temp2, wk, lwork, b_dt, b_dt2, c_a0, local_eqs, eqs, rows,
                          bw, desca, descb, rank, ipiv, info);

        F77NAME(dcopy) (local_eqs, temp, 1, u1, 1);

        // Calculate a(t+1) into temp and copy into a1
        a1_solve_routine(u0, temp, v0, a0, b_dt, b_dt2, c_a0, local_eqs);
        F77NAME(dcopy) (local_eqs, temp, 1, a1, 1);

        // Calculate v(t+1) into temp and copy into v1
        v1_solve_routine(v0, a0, temp, dt, gmm, local_eqs);
        F77NAME(dcopy) (local_eqs, temp, 1, v1, 1);


        // Swap pointers to get f(t+1) into f(t) position
        shift_vec(u0, u1);
        shift_vec(v0, v1);
        shift_vec(a0, a1);

        if (info){
            break;
        }
    }
    // Gather solution onto rank0's F_ref variable
    gather_solution(F_ref, eqs, u0, rank, cores, local_eqs, 0, extra_cols);

    Cblacs_gridexit( ctx );
}

void u1_solve_parallel(double* u0, double* v0, double* a0, double* F, double* M, double* Keff, double* tmp, double* wk,
                       double lwork, double b_dt, double b_dt2, double c_a0, int local_eqs, int eqs, int rows,
                       int bw, int* desca, int* descb, int rank, int* ipiv, int info){

    for (int eq = 0; eq < local_eqs; eq++){
        F[eq] += M[eq] * (u0[eq] * b_dt2 + v0[eq] * b_dt + a0[eq] * c_a0);
    }

    // Create disposable copy of Keff
    F77NAME(dcopy) (local_eqs * rows, Keff, 1, tmp, 1);

    // Solve the system into F
    F77NAME(pdgbsv) (eqs, bw, bw, 1, tmp, 1, desca, ipiv, F, 1, descb, wk, lwork, &info);

    if (info){
        cout << "An error occurred in pdgbsv during PDGBSV: " << info << endl;
    }
}
