int b_addr(int r, int c, int cols);

void print_banded_m(double* M, int rows, int cols);
void print_v(double* M, int rows);

double* clone_vector(double* master, int length);


double* m_diag_add(double* L, double* R, double r_fact, int l, int bw);
// double* vect_vect_multiply(double* L, double* R, int l);
// double* matrix_vector_multiply(double* L, double* R, int l, int bw);
// double* v_v_addition(double* L, double l_fact, double* R, double r_fact int l);
//
