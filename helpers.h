int b_addr(int r, int c, int cols);

void print_banded_m(double* M, int rows, int cols);
void print_v(double* M, int rows);
void print_pos_v(double* M, int cols);

void clone_vector(double* master, double* target, int length);
void shift_vec(double*& pres, double*& futu);
void shift_vec(double*& past, double*& pres, double*& futu);

void m_diag_add(double* L, double* R, double r_fact, int l, int bw);
void multiply_vectors(double* lhs, double* rhs, int eqs);
