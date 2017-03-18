void get_solve_domain(int eqs, int rank, int cores, int begin, int end);
void solve_static(double* K, double* F, int eqs, int bw);
void solve_explicit_parallel(int argc, char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs, int bw);
