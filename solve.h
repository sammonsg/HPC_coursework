void solve_static(double* K, double* F, int eqs, int bw);
void solve_iter(double* u_past, double* u_pres, double* F, double* F_ref, double* KM, double* M, double* M_inv,
                double t, double T_load, int rows, int eqs, int bw);
double solve_explicit(char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs, int bw, double T_load = 1);
double solve_implicit(char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs, int bw, double bt, double gmm,
                    double T_load = 1);
void v1_solve_routine(double* v0, double* a0, double* a1, double dt, double gmm, int eqs);
void a1_solve_routine(double* u0, double* u1, double* v0, double* a0, double b_dt, double b_dt2, double c_a0, int eqs);
