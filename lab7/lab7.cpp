#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#define N 4
#define N_MAX 100

const double TOLX = 1e-8;
const double TOLF = 1e-8;

void wypelnij(double** macierz, double* b, double* x) {
    b[0] = 116.0; b[1] = -226.0; b[2] = 912.0; b[3] = -1174;
    x[0] = 2; x[1] = 2; x[2] = 2; x[3] = 2;
    double zadana_macierz[N][N] = {
        {100.0, -1.0, 2.0, -3.0},
        {1.0, 200.0, -4.0, 5.0},
        {-2.0, 4.0, 300.0, -6.0},
        {3.0, -5.0, 6.0, 400.0}
    };
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            macierz[i][j] = zadana_macierz[i][j];
        }
    }
}

void wyswietl_wektor(double* wektor) {
    for (int i = 0; i < N; i++) {
        cout << setw(20) << wektor[i] << endl;
    }
}

bool koniec(double* residuum, double* estymator_bledu) {
    for (int i = 0; i < N; i++) {
        if (residuum[i] > TOLX || estymator_bledu[i] > TOLF) {
            return false;
        }
    }
    return true;
}

void jacobi(double** macierz, double* b, double* przybl_x) {
    double suma = 0.0;
    double estymator_bledu[N];
    double residuum[N];
    double tmp_x[N];

    printf("NR ======= x =========  RESIDUUM =========== BLAD =====\n");

    for (int k = 0; k < N_MAX; k++) {
        // obliczenie następnego przybliżenia wektora x
        for (int i = 0; i < N; ++i) {
            suma = 0.0;
            for (int j = 0; j < N; ++j) {
                if (i != j) {                    
                    suma += przybl_x[j] * macierz[i][j];
                }    
            }
            tmp_x[i] = (b[i] - suma) / macierz[i][i];
        }

        // obliczenie residuum oraz estymatora bledu
        for (int i = 0; i < N; i++) {
            residuum[i] = 0.0;
            estymator_bledu[i] = fabs(tmp_x[i] - przybl_x[i]);
            // A*x
            for (int j = 0; j < N; j++) {
                residuum[i] += macierz[i][j] * tmp_x[j];
            }
            // A*x - b
            residuum[i] = fabs(residuum[i] - b[i]);
            // nowe przybliżenie
            przybl_x[i] = tmp_x[i];
        }
        
        printf("%1d.=====================================================\n", k+1);
        for (int i = 0; i < N; i++) {
            printf("#  %+.10f |  %.10e | %.10e\n", przybl_x[i], residuum[i], estymator_bledu[i]);
        }
        
        if (koniec(residuum, estymator_bledu)) {
            break;
        }
    }
}

void gauss_seidel(double** macierz, double* b, double* przybl_x) {
    double suma = 0.0;
    double estymator_bledu[N];
    double residuum[N];
    double tmp_x[N];
    double prawa_strona[N];

    printf("NR ======= x =========  RESIDUUM =========== BLAD =====\n");

    for (int k = 0; k < N_MAX; k++) {
        for (int i = 0; i < N; i++) {
            suma = 0.0;
            for (int j = 0; j < N; j++) {
                // U*x
                if (i < j) {
                    suma += macierz[i][j] * przybl_x[j];
                }
            }
            prawa_strona[i] = b[i] - suma;
        } 

        // rozwiazanie rownania 
        for (int i = 0; i < N; i++) {
            suma = 0.0;
            for (int j = 0; j < i; j++) {
                suma += macierz[i][j] * tmp_x[j];
            }
            tmp_x[i] = (prawa_strona[i] - suma) / macierz[i][i];
        }

        for (int i = 0; i < N; i++) {
            residuum[i] = 0.0;
            estymator_bledu[i] = fabs(tmp_x[i] - przybl_x[i]);
            // A*x
            for (int j = 0; j < N; j++) {
                residuum[i] += macierz[i][j] * tmp_x[j];
            }
            // A*x - b
            residuum[i] = fabs(residuum[i] - b[i]);
            // nowe przybliżenie
            przybl_x[i] = tmp_x[i];
        }

        printf("%1d.=====================================================\n", k+1);
        for (int i = 0; i < N; i++) {
            printf("#  %+.10f |  %.10e | %.10e\n", przybl_x[i], residuum[i], estymator_bledu[i]);
        }
        
        if (koniec(residuum, estymator_bledu)) {
            break;
        }        
    }
}

void sor(double** macierz, double* b, double* przybl_x, double omega) {
    double suma;
    double tmp;
    double estymator_bledu[N];
    double residuum[N];
    double tmp_x[N];
    double prawa_strona[N];

    printf("NR ======= x =========  RESIDUUM =========== BLAD =====\n");

    for (int k = 0; k < N_MAX; k++) {
        // obliczenie -[(1-1/omega)D + U]x + b 
        tmp = 1 - 1/omega;
        for (int i = 0; i < N; i++) {
            suma = 0.0;
            for (int j = 0; j < N; j++) {
                if (i <= j) {
                    if (i == j) {
                        suma += tmp * macierz[i][j] * przybl_x[j];
                    } else {
                        suma += macierz[i][j] * przybl_x[j];
                    }
                }
            }
            prawa_strona[i] = b[i] - suma;
        }

        tmp = 1/omega;
        for (int i = 0; i < N; i++) {
            suma = 0.0;
            for (int j = 0; j < i; j++) {
                suma += macierz[i][j] * tmp_x[j];
            }
            tmp_x[i] = (prawa_strona[i] - suma) / (tmp * macierz[i][i]);
        }

        for (int i = 0; i < N; i++) {
            residuum[i] = 0.0;
            estymator_bledu[i] = fabs(tmp_x[i] - przybl_x[i]);
            // A*x
            for (int j = 0; j < N; j++) {
                residuum[i] += macierz[i][j] * tmp_x[j];
            }
            // A*x - b
            residuum[i] = fabs(residuum[i] - b[i]);
            // nowe przybliżenie
            przybl_x[i] = tmp_x[i];
        }

        printf("%1d.=====================================================\n", k+1);
        for (int i = 0; i < N; i++) {
            printf("#  %+.10f |  %.10e | %.10e\n", przybl_x[i], residuum[i], estymator_bledu[i]);
        }

        if (koniec(residuum, estymator_bledu)) {
            break;
        } 
    } 
}

int main() {
    double** A = new double*[N];
    for (int i = 0; i < N; i++) {
        A[i] = new double[N];
    }
    double* b = new double[N];
    double* x = new double[N];

    wypelnij(A, b, x);

    cout << "JACOBI" << endl;
    jacobi(A, b, x);

    cout << endl << endl;

    wypelnij(A, b, x);

    cout << "GAUSS-SEIDEL" << endl;
    gauss_seidel(A, b, x);

    cout << endl << endl;

    wypelnij(A, b, x);

    cout << "SOR" << endl;
    sor(A, b, x, 0.5);

    return 0;
}