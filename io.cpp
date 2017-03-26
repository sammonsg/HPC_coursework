#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

void printHeader(){
    cout <<    "       =============================================================" << endl;
    cout <<    "       =                                                           =" << endl;
    cout <<    "       =  \033[0m HPC Coursework Assignment\033[0m                               =" << endl;
    cout <<    "       =                                                           =" << endl;
    cout <<    "       =  \033[0m \033[1mSubmission by\033[0m: G Sammons                                =" << endl;
    cout <<    "       =  \033[0m \033[1mCID\033[0m: 01252194                                           =" << endl;
    cout <<    "       =                                                           =" << endl;
    cout <<    "       =============================================================" << endl;
}

void printInfo(int N_ele, int eqs, int rank, int cores){
    if (rank==0){
        cout << endl;
        cout << "           " << left << setw(32) << "Number of processes: "         << setw(7) << cores << endl;
        cout << "           " << left << setw(32) << "Global elements:     "         << setw(7) << N_ele << endl;
        cout << "           " << left << setw(32) << "Degrees of freedom:  "         << setw(7) << eqs << endl << endl;
        cout << "       ============================================================="<< endl;
    }
}

void write_v(string filename, double* V, int rows, int cols, int rank = 0, bool zeros = false){
    string path = "output/data/";
    cout << filename << " has been written in " << path << endl;
    std::ofstream stream;
    stream.open (path + filename + ".txt");
    if(zeros) { stream << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl; }
    for (int row = 0; row < rows; row++){
        for (int col = 0; col < cols; col++){
            stream << setprecision(10) << setw(20) << V[row * cols + col];
        }
        stream << endl;
    }
    if(zeros) { stream << setprecision(10) << 0. << setw(20) << 0.  <<  setw(20) << 0.  << endl; }
    stream.close();
}
