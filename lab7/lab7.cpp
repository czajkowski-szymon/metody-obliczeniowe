#include <iostream>

using namespace std;

#define N 4
#define N_MAX 30

const double TOLX = 1e-8;
const double TOLF = 1e-8;

void wyswietl_wektor(double* wektor, int n) {
    for (int i = 0; i < n; i++) {
        cout << setw(20) << wektor[i] << endl;
    }
}

void jacobi(double** macierz, double* b, double* przybl_x) {
    double suma = 0.0;
    double* estymator_bledu = new double[N];
    double* residuum = new double[N];
    double* tmp = new double[N];
    for (int i = 0; i < N_MAX; i++) {
        // obliczenie -(L+U)x + b
        for (int i = 0; i < N; i++) {
            suma = 0.0;
            for (int j = 0; i < N; j++) {
                if (i != j) {
                    suma += przybl_x[j] * macierz[i][j];
                }
            }
            tmp[i] = suma; 
            tmp[i] = b[i] - tmp[i];
        }

        for (int i = 0; i < N; i++) {
            przybl_x[i] = tmp[i] / macierz[i][i];
        }
        
    }
}

void gauss_seidel() {

}

void sor() {

}

int main() {
    double** A = new double*[N];
    for (int i = 0; i < N; i++) {
        A[i] = new double[N];
    }

    double zadana_macierz = {
        {100.0, -1.0, 2.0, -3.0},
        {1.0, 200.0, -4.0, 5.0},
        {-2.0, 4.0, 300.0, -6.0},
        {3.0, -5.0, 6.0, 400.0}
    };

    for (int i = 0; i < N; i++) {
        for (int j = 0; i < N; i++) {
            A[i][j] = zadana_macierz[i][j];
        }
    }

    double* b = new double[N];
    b[0] = 116.0; 
    b[1] = -226.0; 
    b[2] = 912.0;
    b[3] = -1174;

    return 0;
}