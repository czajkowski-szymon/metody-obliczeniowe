#include <iostream>
#include <cmath>
#include <iomanip>

#define N 4 

using namespace std;

void wyswietl_macierz(double** macierz, int n, int* index) {
    for (int i = 0;  i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(10) << setprecision(5) << macierz[index[i]][j];
        }
        cout << endl;
    }
}

void wyswietl_wektor(double* wektor, int n, int* index) {
    for (int i = 0; i < n; i++) {
        cout << setw(5) << wektor[index[i]] << endl;
    }
}

int wybor_elementu_podstawowego(double** macierz, int n, int j, int* index) {
    int nr_wiersza;
    for (int i = j; i < n; i++) {
        if (fabs(macierz[index[i]][j]) < fabs(macierz[index[i+1]][j])) {
            nr_wiersza = index[i+1];
        } else {
            nr_wiersza = index[i];
        }
    }

    return nr_wiersza;
}

/*
 * eliminacja metodą Gaussa
 */
void gauss(double** macierz, int* index, int n) {
    int nr_wiersza;
    double tmp;

    for (int k = 0; k < n - 1; k++) {
        if (macierz[index[k]][k] == 0.0) {
            nr_wiersza = wybor_elementu_podstawowego(macierz, n - 1, index[k], index);
            index[nr_wiersza] = index[k];
            index[k] = nr_wiersza;
        }

        for (int i = k + 1; i < n; i++) {
            tmp = macierz[index[i]][k];
            for (int j = k + 1; j < n; j++) {
                macierz[index[i]][j] = macierz[index[i]][j] - macierz[index[k]][j] * (tmp / macierz[index[k]][k]);
            }  
            macierz[index[i]][k] = tmp / macierz[index[k]][k];
        }
    }
}

/*
 * rozwiązanie równania Lx = y
 */
void oblicz_y(double** macierz, double* wektor, int n, int* index) {
    double suma = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            suma += macierz[index[i]][j] * wektor[index[j]];
        }
        wektor[index[i]] = (wektor[index[i]] - suma) / 1.0;
        suma = 0.0;
    }
}

/*
 * rozwiązanie równania Ux = b
 */
void oblicz_x(double** macierz, double* wektor, int n, int* index) {
    double suma = 0.0;
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            suma += macierz[index[i]][j] * wektor[index[j]];
        }
        wektor[index[i]] = (wektor[index[i]] - suma) / macierz[index[i]][i];
        suma = 0.0;
    }
}


int main() {
    // tworzenie macierzy A
    double** A = new double*[N];
    for (int i = 0; i < N; i++) {
        A[i] = new double[N];
    }

    double zadana_macierz[N][N] = {
        {1.0, -20.0, 30.0, -4.0},
        {2.0, -40.0, -6.0, 50.0},
        {9.0, -180.0, 11.0, -12.0},
        {-16.0, 15.0, -140.0, 13.0}
    };

    // wypełnienie macierzy A danymi
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = zadana_macierz[i][j];
        }
    }

    // stworzenie wektora b i wypełnienie go danymi
    double* b = new double[N];
    b[0] = 35.0; 
    b[1] = 104.0;
    b[2] = -366.0;
    b[3] = -354.0;

    // wektor indeksów
    int index[4] = {0, 1, 2, 3};

    gauss(A, index, N);
    cout << "Polaczona macierz L i U:" << endl;
    wyswietl_macierz(A, N, index);
    oblicz_y(A, b, N, index);
    oblicz_x(A, b, N, index);
    cout << "Wektor wyniku:" << endl;
    wyswietl_wektor(b, N, index);

    for (int i = 0; i < N; i++) {
        delete[] A[i];
    }
    delete[] A;

    return 0;
}