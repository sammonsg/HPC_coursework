using namespace std;

void mk_k_mat(double* K_ele, double* K, int n, int dof, int rows);
void mk_k_ele(char* argv[], double l, double* K, int rows);
void mk_k_ele1(char* argv[], double l, double* K, int rows);
int addr(int r, int c, int cols);

void mk_banded_k_mat(double* K_ele, double* K, int n, int dof, int rows, int bw);
void mk_banded_k_ele(char* argv[], double l, double* K, int rows, int bw);
void mk_banded_k_ele1(char* argv[], double l, double* K, int rows);
int b_addr(int r, int c, int bw);


void mk_F_mat(double* F_ele, double* F, int n, int dof, int rows);
void mk_F_ele(double* F, double q_x, double q_y, double l);
