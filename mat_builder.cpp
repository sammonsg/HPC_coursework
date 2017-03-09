#include <cmath>
#include <cstdlib>
#include <iostream>
#include "helpers.h"

using namespace std;


void mk_banded_k_mat(double* K_ele, double* K, int N, int dof, int rows, int bw){
    int K_ele_len = rows * 3;
    // cout << K_ele_len << " and " << N << endl;
    for (int node = 0; node < N; node++){
        for (int element = 0; element < K_ele_len; element++){
            // cout << element + node * K_ele_len << endl;
            K[element + node * K_ele_len] = K_ele[element];
        }
    }
}

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

void mk_F_ele(double* F, double q_x, double q_y, double l, int dof){
    double * F_ele = new double[dof*2]();

    // From a computational point of view, this step is unnecessary due to both
    // boundary condition fully restrained, but if this was not the case, the
    // elemental forces are available as below
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

void mk_F_mat(double* F, int N, int dof, double q_x, double q_y, double l){
    double * F_node = new double[dof]();
    mk_F_ele(F_node, q_x, q_y, l, dof);
    for (int a = 0; a < N; a++){
        for (int dir = 0; dir < dof; dir++){
                F[dir + a * dof] += F_node[dir];
        }
    }
    delete [] F_node;
}
