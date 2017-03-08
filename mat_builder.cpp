#include <cmath>
#include <cstdlib>
#include <iostream>
#include "helpers.h"

using namespace std;

// void mk_k_mat(double* K_ele, double* K, int n, int dof, int rows){
//     for (int a = -1; a < n + 1; a++){
//         for (int r = 1; r < dof*2+1; r++){
//             for (int c = 1; c < dof*2+1; c++){
//                 int col = r + a * dof;
//                 int row = c + a * dof;
//                 if(min(row, col) > -1 && max(row,col) < rows+1){
//                     K[addr(row,col,rows)] += K_ele[addr(r,c,dof*2)];
//                 }
//             }
//         }
//     }
//     print_m(K, rows);
// }
//
// void mk_k_ele(char* argv[], double l, double* K, int rows){
//     // Units are converted from those supplied to SI units
//     double A = atof(argv[3]) * pow(10, -6);
//     double I = atof(argv[4]) * pow(10, -12);
//     double E = atof(argv[5]) * pow(10, 6);
//     // K will be stored in lower-diagonal form, so will be the elemental form
//     K[addr(1,1,rows)] = A*E/l;
//     K[addr(4,1,rows)] = -K[addr(1,1,rows)];
//     K[addr(2,2,rows)] = 12 * E * I / pow(l, 3);
//     K[addr(3,2,rows)] = 6 * E * I / pow(l, 2);
//     K[addr(5,2,rows)] = -K[addr(2,2,rows)];
//     K[addr(6,2,rows)] = K[addr(3,2,rows)];
//     K[addr(3,3,rows)] = 4 * E * I / l;
//     K[addr(5,3,rows)] = -K[addr(3,2,rows)];
//     K[addr(6,3,rows)] = 2 * E * I / l;
//     K[addr(4,4,rows)] = K[addr(1,1,rows)];
//     K[addr(5,5,rows)] = K[addr(2,2,rows)];
//     K[addr(6,5,rows)] = K[addr(5,3,rows)];
//     K[addr(6,6,rows)] = K[addr(3,3,rows)];
//     print_m(K, 6);
// }
//
// void mk_k_ele1(char* argv[], double l, double* K, int rows){
//     // Make ones everywhere
//     for (int r = 1; r < rows+1; r++){
//         for (int c = 1; c < rows+1; c++){
//             K[addr(r,c,rows)] = 1;
//         }
//     }
// }

void mk_banded_k_mat(double* K_ele, double* K, int N, int dof, int rows, int bw){
    int K_ele_len = rows * 3;
    // cout << K_ele_len << " and " << N << endl;
    for (int node = 0; node < N; node++){
        for (int element = 0; element < K_ele_len; element++){
            // cout << element + node * K_ele_len << endl;
            K[element + node * K_ele_len] = K_ele[element];
        }
    }
    // print_m(K, rows);
}

void mk_banded_k_ele(char* argv[], double l, double* K_ele, int rows, int bw){
    // Make temporary 13 x 6 matrix
    double * K = new double[rows * 3];

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
}

void mk_F_mat(double* F_ele, double* F, int n, int dof, int rows){
    for (int a = 0; a < n; a++){
        for (int r = 1; r < dof*2+1; r++){
                // F[r + a * dof)] += F_ele[addr(r)];
        }
    }
}

void mk_F_ele(double* F, double q_x, double q_y, double l){
    // Units are converted from those supplied to SI units
    F[0] = 0.5 * q_x * l;
    F[1] = 0.5 * q_y * l;
    F[2] = q_y * l * l / 12.0;
    F[3] = F[0];
    F[4] = F[1];
    F[5] = -F[2];
}
