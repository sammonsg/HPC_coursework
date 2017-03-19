void mk_banded_k_mat(char* argv[], double l, double* K, int N, int dof, int rows, int bw);
void mk_banded_k_ele(char* argv[], double l, double* K, int rows, int bw);
void mk_M_mat(char* argv[], double* M, int N, double l, int dof);

void mk_F_mat(double* F, int N, int dof, double q_x, double q_y, int F_centre,  double l);
void mk_F_ele(double* F, double q_x, double q_y, double l, int dof);

void mk_km_mat(double* KM, double* K_ref, double* M, int eqs, int bw);
void mk_keff_mat(double* KM, double* K_ref, double* M, double b_dt2, int eqs, int bw);

void mk_truncated_mat(double* K, double* K_ref, int rows, int begin, int end, int row_offset);

void invert_v(double* out, double* in, int n);
