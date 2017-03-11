int b_addr(int r, int c, int cols);

void print_banded_m(double* M, int rows, int cols);
void print_v(double* M, int rows);

double* clone_vector(double* master, int length);

void vect_vect_multiply(double* L, double* R, double* out, int l);
void banded_vector_multiply(double* L, double* R, double* out, int l, int bw);
void vect_vect_addition(double* L, double* R, double* out, int l, int sign);
