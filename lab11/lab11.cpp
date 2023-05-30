#include <iostream>
#include <fstream>
#include <cmath>
#include "calerf.h"

using namespace std;

// parametry
const double T_MIN = 0.0;
const double T_MAX = 2.0;
const double D = 1.0;

const double POCZATEK_PRZEDZIALU = -9.0;
const double KONIEC_PRZEDZIALU = 9.0;

//const double h = 0.1;

const double P_LAMBDA = 1.0;

double** stworz_macierz(int N, int M) {
    double** macierz = new double*[N];
    for (int i = 0; i < N; i++) {
        macierz[i] = new double[M];
    }
    return macierz;
}

void usun_macierz(double** macierz, int N) {
    for (int i = 0; i < N; i++) {
        delete[] macierz[i];
    }
    delete[] macierz;
}


double* siatka_kroku(double h, int M) {
    double* wektor = new double[M];
    double x = POCZATEK_PRZEDZIALU;

    for (int i = 0; i < M; i++) {
        wektor[i] = x;
        x += h;
    }

    return wektor;
}

void zapisz_do_pliku_w(double* wektor, int N, string nazwa) {
    ofstream plik;
    plik.open(nazwa);
    for (int i = 0; i < N; i++) {
        plik << wektor[i] << endl;
    }   
}


void zapisz_do_pliku_m(double** macierz, int N, int M, string nazwa) {
    ofstream plik;
    plik.open(nazwa);
    double* kroki = siatka_kroku(0.1, N);
    for (int i = 0; i < N; i++) {
        plik << kroki[i] << "\t";
        for (int j = 0; j < M; j++) {
            plik << macierz[i][j] << "\t";
        }
        plik << endl;
    }

    plik.close();
}

double rozwiazanie_analityczne(double x, double t) {
    return 0.5 * calerfpack::erfc_l(x / 2 * sqrt(D * t));
}

double** oblicz_rozwiazanie_analityczne(double h, double dt, int N, int M) {
    double** macierz_rozwiazan = stworz_macierz(N, M);
    double x = POCZATEK_PRZEDZIALU;
    double t = T_MIN;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            macierz_rozwiazan[i][j] = rozwiazanie_analityczne(x, t);
            x += h;
        }
        x = POCZATEK_PRZEDZIALU;
        t += dt;
    }

    return macierz_rozwiazan;
}

void warunek_poczatkowy(double** macierz, double h, int N, int M) {
    double x = POCZATEK_PRZEDZIALU;
    for (int i = 0; i < M; i++) {
        if (x < 0.0) {
            macierz[0][i] = 1.0;
        } else {
            macierz[0][i] = 0.0;
        }
        x += h;
    }
}

void warunek_brzegowy(double** macierz, int N, int M) {
    for (int i = 0; i < N; i++) {
        macierz[i][0] = 1.0;
        macierz[i][M - 1] = 0.0;
    }
}

void thomas(double *l, double *d, double *u, double *b, double *x, int M) {
    for (int i = 2; i < M; i++) {
        d[i] = d[i] - (l[i - 1] / d[i - 1]) * u[i - 1];
        b[i] = b[i] - (l[i - 1] / d[i - 1]) * b[i - 1];
    }

    x[M - 1] = b[M - 1] / d[M - 1];

    for (int i = M - 2; i >= 0; i--) {
        x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
    }
}

// void gauss_seidel(double** macierz, double* b, double* przybl_x) {
//     double suma = 0.0;
//     double estymator_bledu[N];
//     double residuum[N];
//     double tmp_x[N];
//     double prawa_strona[N];

//     for (int k = 0; k < N_MAX; k++) {
//         for (int i = 0; i < N; i++) {
//             suma = 0.0;
//             for (int j = 0; j < N; j++) {
//                 // U*x
//                 if (i < j) {
//                     suma += macierz[i][j] * przybl_x[j];
//                 }
//             }
//             prawa_strona[i] = b[i] - suma;
//         } 

//         // rozwiazanie rownania 
//         for (int i = 0; i < N; i++) {
//             suma = 0.0;
//             for (int j = 0; j < i; j++) {
//                 suma += macierz[i][j] * tmp_x[j];
//             }
//             tmp_x[i] = (prawa_strona[i] - suma) / macierz[i][i];
//         }

//         for (int i = 0; i < N; i++) {
//             residuum[i] = 0.0;
//             estymator_bledu[i] = fabs(tmp_x[i] - przybl_x[i]);
//             // A*x
//             for (int j = 0; j < N; j++) {
//                 residuum[i] += macierz[i][j] * tmp_x[j];
//             }
//             // A*x - b
//             residuum[i] = fabs(residuum[i] - b[i]);
//             // nowe przybliżenie
//             przybl_x[i] = tmp_x[i];
//         }
        
//         if (koniec(residuum, estymator_bledu)) {
//             break;
//         }        
//     }
// }

void laasonen_thomas_pom(double** macierz, int N, int M) {
    double lambda = 1.0 + 2.0*P_LAMBDA;
    double* u = new double[M];
    double* d = new double[M];
    double* l = new double[M];
    double* b = new double[M];
    double* x = new double[M];

    for (int i = 1; i < N; i++) {
        u[0] = 0.0;
        d[0] = 1.0;
        l[0] = 0.0;
        b[0] = macierz[i - 1][0];
        for (int j = 1; j < M - 1; j++) {
            u[j] = P_LAMBDA;
            d[j] = -lambda;
            l[j] = P_LAMBDA;
            b[j] = -macierz[i - 1][j];
        }
        u[M - 1] = 0.0;
        d[M - 1] = 1.0;
        l[M - 1] = 0.0;
        b[M - 1] = 0.0;

        thomas(l, d, u, b, x, M);

        for (int j = 1; j < M - 1; j++) {
            macierz[i][j] = x[j];
        }
    }

    delete[] u;
    delete[] d;
    delete[] l;
    delete[] b;
    delete[] x;
}

// void laasonen_gs_pom(double** macierz, int N, int M) {
//     double lambda = 1.0 + 2.0*P_LAMBDA;
//     double* b = new double[m];
//     double* wyn = new double[m];
    
//     for (int i = 0; i < m; i++) {
//         b[i] = 0 .0;
//         wyn[i] = 0 .0;
//     }

//     double** new_A = new double*[m];
//     for ( i nt i = 0 ; i < m ; i ++) {
//         new_A[i] = new double[m];
//     }

//     for (int i = 0; i < m; i++) {
//         for (int j = 0; j < m; j++) {
//             new_A[i][j] = 0 .0;
//         }
//     }
//     for (int k = 1; k < n; k++) {
//         new_A[0][0] = 1 .0;
//         b[0] = A[k - 1][0];
//         for (int i = 1; i < m - 1; i++) {
//             new_A[i][i] = - lambda; 
//             new_A[i][i + 1] = P_LAMBDA; 
//             new_A[i][i - 1] = P_LAMBDA; 
//             b[i] = - A[k - 1][i];
//         }

//         b[m - 1] = 0.0;
//         new_A[m - 1][m - 1] = 1.0;
//         gauss_seidel(new_A, b, wyn, n, m);
//         for (int i = 1; i < m - 1; i++) {
//             A[k][i] = wyn[i];
//         }
//     }
// }

double** policz_bledy(double** wyniki_analityczne, double** wyniki_numeryczne, int N, int M) {
    double** bledy = stworz_macierz(N, M);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            bledy[i][j] = fabs(wyniki_numeryczne[i][j] - wyniki_analityczne[i][j]);
        }
    }

    return bledy;
}

double norma_maksimum(double* wektor, int N) {
    double max = fabs(wektor[0]);
    for (int i = 1; i < N; i++) {
        if (max < fabs(wektor[i])) {
            max = wektor[i];
        }
    }
    return max;
}

double* bledy_max(double** bledy, int N, int M) {
    double* bledy_max = new double[N];
    for (int i = 0; i < N; i++) {
        bledy_max[i] = norma_maksimum(bledy[i], M);
    }
    return bledy_max;
}

double** laasonen_thomas(double h, int N, int M) {
    double** laasonen = stworz_macierz(N, M);
    warunek_poczatkowy(laasonen, h, N, M);
    warunek_brzegowy(laasonen, N, M);
    laasonen_thomas_pom(laasonen, N, M);
    return laasonen;
}

// double** laasonen_gs(double h, int N, int M) {
//     double** laasonen = stworz_macierz(N, M);
//     warunek_poczatkowy(laasonen, h, N, M);
//     warunek_brzegowy(laasonen, N, M);
//     laasonen_gs_pom(laasonen, N, M);
//     return laasonen;
// }

double** transponuj(double** macierz, int N, int M) {
    double** result = stworz_macierz(M, N);
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            result[j][i] = macierz[i][j];
        }
    }

    return result;
}

int main() {
    /*
     *  LAASONEN - THOMAS 
     */
    double h = 0.5;
    double dt = (P_LAMBDA*h*h) / D;
    int N = (T_MAX - T_MIN) / dt;
    int M = (KONIEC_PRZEDZIALU - POCZATEK_PRZEDZIALU) / h;
    
    double** rozwiazanie_analityczne;
    // zapisz_do_pliku_m(transponuj(rozwiazanie_analityczne, N, M), M, N, "rozwiazanie_analityczne.txt");
    double** laasonen_thomas_sol;
    // zapisz_do_pliku_m(transponuj(laasonen_thomas_sol, N, M), M, N, "laasonen_thomas_numerycznie.txt");
    double** bledy;
    double maks_blad;

    ofstream bledy_plik;
    bledy_plik.open("lt_bledy.txt");
    for (h; h > 1e-2; h /= 2) {
        double dt = (P_LAMBDA*h*h) / D;
        int N = (T_MAX - T_MIN) / dt;
        int M = (KONIEC_PRZEDZIALU - POCZATEK_PRZEDZIALU) / h;
        rozwiazanie_analityczne = oblicz_rozwiazanie_analityczne(h, dt, N, M);
        laasonen_thomas_sol = laasonen_thomas(h, N, M);
        bledy = policz_bledy(rozwiazanie_analityczne, laasonen_thomas_sol, N, M);
        maks_blad = norma_maksimum(bledy[N-1], M);
        bledy_plik << h << "\t" << maks_blad << endl;
    }

    return 0;
}