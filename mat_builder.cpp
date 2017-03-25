#include <cmath>
#include <cstdlib>
#include <iostream>
#include "helpers.h"

using namespace std;

#include "lapack_c.h"
#include "helpers.h"


void mk_banded_k_ele(char* argv[], double l, double* K_ele, int rows, int bw){
    // Make temporary 13 x 6 matrix
    double * K = new double[rows * 6]();

    // Units are converted from those supplied to SI units
    double A = atof(argv[3]) * pow(10, -6);
    double I = atof(argv[4]) * pow(10, -12);
    double E = atof(argv[5]) * pow(10, 6);


    // K will be stored in lower-diagonal form, so will be the elemental form
    K[b_addr(1,1,bw)] = A*E/l;
    K[b_addr(4,1,bw)] = -K[b_addr(1,1,bw)];
    K[b_addr(2,2,bw)] = 12 * E * I / pow(l, 3);
    K[b_addr(3,2,bw)] = 6 * E * I / pow(l, 2);
    K[b_addr(5,2,bw)] = -K[b_addr(2,2,bw)];
    K[b_addr(6,2,bw)] = K[b_addr(3,2,bw)];
    K[b_addr(3,3,bw)] = 4 * E * I / l;
    K[b_addr(5,3,bw)] = -K[b_addr(3,2,bw)];
    K[b_addr(6,3,bw)] = 2 * E * I / l;
    K[b_addr(4,4,bw)] = K[b_addr(1,1,bw)];
    K[b_addr(5,5,bw)] = K[b_addr(2,2,bw)];
    K[b_addr(6,5,bw)] = K[b_addr(5,3,bw)];
    K[b_addr(6,6,bw)] = K[b_addr(3,3,bw)];

    // Add symmetric terms
    for (int c = 2; c < 7; c++){
        for (int r = 1; r < c; r++){
            K[b_addr(r,c,bw)] = K[b_addr(c,r,bw)];
        }
    }

    // Add columns 1,2,3 to 4,5,6 to compose non-edge node equation block
    int K_ele_len = (1 + 3 * bw) * 3;
    for (int n = 0; n < K_ele_len; n++){
        K_ele[n] = K[n] + K[n + K_ele_len];
    }
    delete [] K;
}


void mk_banded_k_mat(char* argv[], double l, double* K, int N, int dof, int rows, int bw){
    int K_ele_len = rows * dof;
    double * K_ele = new double[K_ele_len]();
    mk_banded_k_ele(argv, l, K_ele, rows, bw);

    // cout << K_ele_len << " and " << N << endl;
    for (int node = 0; node < N; node++){
        for (int element = 0; element < K_ele_len; element++){
            // cout << element + node * K_ele_len << endl;
            K[element + node * K_ele_len] = K_ele[element];
        }
    }
    delete [] K_ele;
}


void mk_F_ele(double* F, double q_x, double q_y, double l, int dof){
    double * F_ele = new double[dof*2]();

    // From a computational point of view, this step is unnecessary due to both boundary condition fully restrained,
    // but if this was not the case, the elemental forces are available as below
    F_ele[0] = 0.5 * q_x * l;
    F_ele[1] = 0.5 * q_y * l;
    F_ele[2] = q_y * l * l / 12.0;
    F_ele[3] = F_ele[0];
    F_ele[4] = F_ele[1];
    F_ele[5] = -F_ele[2];

    // As all free nodes have a left and right element, add nodal forces
    for (int a = 0; a < dof; a++){
        F[a] = F_ele[a] + F_ele[a + dof];
    }
    delete [] F_ele;

}

void mk_F_mat(double* F, int N, int dof, double q_x, double q_y, int F_centre, double l){
    double * F_node = new double[dof]();
    mk_F_ele(F_node, q_x, q_y, l, dof);
    for (int a = 0; a < N; a++){
        for (int dir = 0; dir < dof; dir++){
                F[dir + a * dof] += F_node[dir];
        }
    }
    F[(N * dof - 1) / 2] += F_centre;

    delete [] F_node;
}

void mk_M_mat(char* argv[], double* M, int N, double l, int dof){
    for (int n = 0; n < N; n++){
        // Inertia terms are doubled due to overlay
        M[0 + n * dof] = 1;
        M[1 + n * dof] = 1;
        M[2 + n * dof] = l * l / 12.0;
    }
}

void mk_km_mat(double* KM, double* K_ref, double* M, int eqs, int bw){
    int banded_rows = 3 * bw + 1;
    int rows = 2 * bw + 1;
    double * KM_full = new double[eqs * banded_rows];
    F77NAME(dcopy) (banded_rows * eqs, K_ref, 1, KM_full, 1);
    for (int c = 0; c < eqs; c++){
        for (int r = 0; r < rows; r++){
            KM[r + c * rows] = KM_full[r + 4 + c * banded_rows];
        }
        KM[c * rows + bw] += -2 * M[c];
    }
    delete [] KM_full;
}

void mk_keff_mat(double* Keff, double* K_ref, double* M, double b_dt2, int eqs, int bw, int offset = 0){
    int rows = 3 * bw + 1 + offset;
    F77NAME(dcopy) (rows * eqs, K_ref, 1, Keff, 1);
    for (int c = 0; c < eqs; c++){
        Keff[2 * bw + offset + c * rows] += M[c] * b_dt2;
    }
}


void mk_truncated_mat(double* K, double* K_ref, int rows, int begin, int end, int row_offset){
    int ref_rows = rows + row_offset;

    for (int c = begin; c < end + 1; c++){
        for (int r = 0; r <  rows - max(-row_offset,0); r++){
            K[r + max(-row_offset,0) + (c-begin) * rows] = K_ref[r + max(row_offset,0) + c * ref_rows];
        }
    }
}

void mk_truncated_v(double* K, double* K_ref, int begin, int end){
    int len = end - begin + 1;
    for (int r = 0; r < len; r++){
        K[r] = K_ref[r + begin];
    }
}

void invert_v(double* out, double* in, int n){
    for (int a = 0; a < n; a++){
        out[a] = 1.0 / in[a];
    }
}
