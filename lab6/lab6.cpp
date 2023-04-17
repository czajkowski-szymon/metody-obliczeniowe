#include <iostream>
#include <iomanip>

using namespace std;

#define N 6

void oblicz_a(double* eta, double* u, double* d, double* l) {
    // pierwszy krok
    eta[0] = d[0];

    // dla reszty
    for (int i = 1; i < N; i++) {
        eta[i] = d[i] - (l[i - 1] * u[i - 1]) / eta[i - 1];
    }
}

void oblicz_b(double* eta, double* r, double* b, double* l) {
    // pierwszy krok
    r[0] = b[0];

    // dla reszty
    for (int i = 1; i < N; i++) {
        r[i] = b[i] - (l[i - 1] * r[i - 1]) / eta[i - 1];
    }
}

void oblicz_x(double* x, double* eta, double* r, double* u) {
    // pierwszy krok
    x[N-1] = r[N-1] / eta[N-1];

    // dla reszty
    for (int i = N - 2; i >= 0; i--) {
        x[i] = (r[i] - u[i] * x[i + 1]) / eta[i];
    }
}

int main() {
    double* u = new double[N - 1];
    double* d = new double[N];
    double* l = new double[N - 1];
    double* b = new double[N];
    double* eta = new double[N];
    double* r = new double[N];
    double* x = new double[N];

    u[0] = 1.0 / 2.0; 
    u[1] = 1.0 / 4.0;
    u[2] = 1.0 / 6.0;
    u[3] = 1.0 / 8.0;
    u[4] = 1.0 / 10.0;

    d[0] = 10.0;
    d[1] = 20.0;
    d[2] = 30.0;
    d[3] = 30.0;
    d[4] = 20.0;
    d[5] = 10.0;
    
    l[0] = 1.0 / 3.0;
    l[1] = 1.0 / 5.0;
    l[2] = 1.0 / 7.0;
    l[3] = 1.0 / 9.0;
    l[4] = 1.0 / 11.0;

    b[0] = 31.0;
    b[1] = 165.0 / 4.0;
    b[2] = 917.0 / 30.0;
    b[3] = 851.0 / 28.0;
    b[4] = 3637.0 / 90.0;
    b[5] = 332.0 / 11.0;

    oblicz_a(eta, u, d, l);
    oblicz_b(eta, r, b, l);
    oblicz_x(x, eta, r, u);

    for (int i = 0; i < N; i++) {
        cout << "x" << i + 1 << " = " << setprecision(10) <<  x[i] << endl;
    }

    delete[] u;
    delete[] d;
    delete[] l;
    delete[] b;
    delete[] eta;
    delete[] r;

    return 0;
}