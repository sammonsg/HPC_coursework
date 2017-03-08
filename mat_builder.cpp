#include <cmath>
#include <cstdlib>
#include <iostream>
#include "helpers.h"

using namespace std;

void mk_k_mat(double* K_ele, double* K, int n, int dof, int rows){
    for (int a = -1; a < n + 1; a++){
        for (int r = 1; r < dof*2+1; r++){
            for (int c = 1; c < dof*2+1; c++){
                int col = r + a * dof;
                int row = c + a * dof;
                if(min(row, col) > -1 && max(row,col) < rows+1){
                    K[addr(row,col,rows)] += K_ele[addr(r,c,dof*2)];
                }
            }
        }
    }
    print_m(K, rows);
}

void mk_k_ele(char* argv[], double l, double* K, int rows){
    // Units are converted from those supplied to SI units
    double A = atof(argv[3]) * pow(10, -6);
    double I = atof(argv[4]) * pow(10, -12);
    double E = atof(argv[5]) * pow(10, 6);
    // K will be stored in lower-diagonal form, so will be the elemental form
    K[addr(1,1,rows)] = A*E/l;
    K[addr(4,1,rows)] = -K[addr(1,1,rows)];
    K[addr(2,2,rows)] = 12 * E * I / pow(l, 3);
    K[addr(3,2,rows)] = 6 * E * I / pow(l, 2);
    K[addr(5,2,rows)] = -K[addr(2,2,rows)];
    K[addr(6,2,rows)] = K[addr(3,2,rows)];
    K[addr(3,3,rows)] = 4 * E * I / l;
    K[addr(5,3,rows)] = -K[addr(3,2,rows)];
    K[addr(6,3,rows)] = 2 * E * I / l;
    K[addr(4,4,rows)] = K[addr(1,1,rows)];
    K[addr(5,5,rows)] = K[addr(2,2,rows)];
    K[addr(6,5,rows)] = K[addr(5,3,rows)];
    K[addr(6,6,rows)] = K[addr(3,3,rows)];
    print_m(K, 6);
}

void mk_k_ele1(char* argv[], double l, double* K, int rows){
    // Make ones everywhere
    for (int r = 1; r < rows+1; r++){
        for (int c = 1; c < rows+1; c++){
            K[addr(r,c,rows)] = 1;
        }
    }
}

void mk_banded_k_mat(double* K_ele, double* K, int n, int dof, int rows, int bw){
    for (int a = -1; a < n + 1; a++){
        for (int r = 1; r < dof*2+1; r++){
            for (int c = 1; c < dof*2+1; c++){
                int col = r + a * dof;
                int row = c + a * dof;
                if(min(row, col) > -1 && max(row,col) < rows+1){
                    K[addr(row,col,rows)] += K_ele[addr(r,c,dof*2)];
                }
            }
        }
    }
    print_m(K, rows);
}

void mk_banded_k_ele(char* argv[], double l, double* K, int rows, int bw){
    // Units are converted from those supplied to SI units
    double A = atof(argv[3]) * pow(10, -6);
    double I = atof(argv[4]) * pow(10, -12);
    double E = atof(argv[5]) * pow(10, 6);
    // K will be stored in lower-diagonal form, so will be the elemental form
    K[b_addr(1,1,rows)] = A*E/l;
    K[b_addr(4,1,rows)] = -K[b_addr(1,1,rows)];
    K[b_addr(2,2,rows)] = 12 * E * I / pow(l, 3);
    K[b_addr(3,2,rows)] = 6 * E * I / pow(l, 2);
    K[b_addr(5,2,rows)] = -K[b_addr(2,2,rows)];
    K[b_addr(6,2,rows)] = K[b_addr(3,2,rows)];
    K[b_addr(3,3,rows)] = 4 * E * I / l;
    K[b_addr(5,3,rows)] = -K[b_addr(3,2,rows)];
    K[b_addr(6,3,rows)] = 2 * E * I / l;
    K[b_addr(4,4,rows)] = K[b_addr(1,1,rows)];
    K[b_addr(5,5,rows)] = K[b_addr(2,2,rows)];
    K[b_addr(6,5,rows)] = K[b_addr(5,3,rows)];
    K[b_addr(6,6,rows)] = K[b_addr(3,3,rows)];
    for (int r = 1; r < 6)
}

void mk_banded_k_ele1(char* argv[], double l, double* K, int rows, int bw){
    // Make ones everywhere
    for (int r = 1; r < rows+1; r++){
        for (int c = 1; c < rows+1; c++){
            K[addr(r,c,rows)] = 1;
        }
    }
}

void mk_F_mat(double* F_ele, double* F, int n, int dof, int rows){
    for (int a = 0; a < n; a++){
        for (int r = 1; r < dof*2+1; r++){
                F[addr(r + a * dof)] += F_ele[addr(r)];
        }
    }
}

void mk_F_ele(double* F, double q_x, double q_y, double l){
    // Units are converted from those supplied to SI units
    F[addr(1)] = 0.5 * q_x * l;
    F[addr(2)] = 0.5 * q_y * l;
    F[addr(3)] = q_y * l * l / 12.0;
    F[addr(4)] = F[addr(1)];
    F[addr(5)] = F[addr(2)];
    F[addr(6)] = -F[addr(3)];
}
