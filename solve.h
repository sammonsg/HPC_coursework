void solve_static(double* K, double* F, int eqs, int bw);
void solve_iter(double* u_past, double* u_pres, double* F, double* F_ref, double* KM, double* M, double* M_inv,
                double t, int rows, int eqs, int bw);
void solve_explicit(char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs, int bw);
void solve_implicit(char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs, int bw, double bt, double gmm);
