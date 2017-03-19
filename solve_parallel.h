void get_solve_domain(int eqs, int rank, int cores, int* begin, int* end);
void solve_explicit_parallel(char* argv[], double* K_ref, double* F_ref, double* M_ref, int eqs, int bw,
                            int rank, int cores);
