#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

const double P = 1.0, Q = 0.0, R = -4.0;

const double ALFA = 0.0, BETA = 1.0, GAMMA = -1.0;
const double FI = 0.0, PSI = 1.0, TETA = 0.0;

const double A = 0.0, B = 1.0;

double rozwiazanie_analityczne(double x) {
    return (exp(2.0 - 2.0 * x) - 4.0 * exp(4.0 - x * 2.0) + 4.0 * exp(x * 2.0) - exp(2.0 + 2.0 * x) - x + x * exp(4.0))
        / (4.0 - 4 * exp(4.0));
}

void thomas(double *l, double *d, double *u, double *b, double *x, int N) {
    double *r = new double[N];
    double *eta = new double[N];

    // Realizuje algorytm Thomasa zgodnie z wykladem
    eta[0] = d[0];
    r[0] = b[0];

    // wyliczenie eta
    for (int i = 1; i < N; i++) {
        eta[i] = d[i] - l[i - 1] * (u[i - 1] / eta[i - 1]);
    }

    // wyliczenie r
    for (int i = 1; i < N; i++) {
        r[i] = b[i] - l[i - 1] * r[i - 1] / eta[i - 1];
    }

    // Obliczanie rozwiazania
    x[N - 1] = r[N - 1] / eta[N - 1];

    // Obliczenie pozostalych x
    for (int i = N - 2; i >= 0; i--) {
        x[i] = (r[i] - u[i] * x[i + 1]) / eta[i];
    }

    delete[] r;
    delete[] eta;
}

int norma_maksimum(double *blad, int N) {
    double maksymalny = fabs(blad[0]);    // Inicjalizujemy blad pierwsza wartoscia wektora
    int index = 0;

    for (int i = 0; i < N; i++)
        if (fabs(blad[i]) > maksymalny) {
        maksymalny = fabs(blad[i]);
        index = i;
        }

    return index;
}

double numerow(double h, int N) {
    double *l, *d, *u, *b, *x, *blad, x0 = A, xn = A;
    fstream numerow, analitycznie;
    numerow.open("wyniki_numerow.txt", ios_base::app);
    analitycznie.open("wyniki_analityczne.txt", ios_base::app);
    analitycznie << scientific;
    numerow << scientific;
    cout.precision(10);
    l = new double[N];
    d = new double[N];
    u = new double[N];
    b = new double[N];
    x = new double[N];
    blad = new double[N];

    // Realizuje algorytm zgodnie ze wzorami z wykladu

    // Z warunku brzegowego
    u[0] = ALFA / h;
    d[0] = BETA - ALFA / h;
    b[0] = -GAMMA;

    // Wyznaczenie srodkowych wyrazow w zapisie wektorowym
    for (int i = 1; i < N - 1; i++) {
        l[i - 1] = P / (h * h) + R / 12.0;
        d[i] = (-2.0 * P) / (h * h) + R * (10.0 / 12.0);
        u[i] = P / (h * h) + R / 12.0;
        b[i] = (x0 + i * h - h) / 12.0 + (10.0 / 12.0) * (x0 + i * h) + (x0 + i * h + h) / 12.0;
    }

    // Z warunku brzegowego
    l[N - 2] = -FI / h;
    d[N - 1] = -FI / h + PSI;
    b[N - 1] = -TETA;


    // Rozwiazanie macierzy trojdiagonalnej algorytmem thomasa
    thomas(l, d, u, b, x, N);

    // Obliczenie bledu między algotytmem numerow, a rozwiązaniem analitycznym
    for (int i = 0; i < N; i++) {
        blad[i] = fabs(x[i] - rozwiazanie_analityczne(xn));
        xn += h;
    }

    int naj = norma_maksimum
(blad, N); //znajdujemy największy błąd

    if (N == 1002) {
        for (int i = 0; i < N; i++) {
        numerow << x0 << "\t" << x[i] << "\t\n";
        analitycznie << x0 << "\t" << rozwiazanie_analityczne(x0) << "\n";
        x0 += h;
        }
    }

    delete[] l;
    delete[] d;
    delete[] u;
    delete[] x;
    delete[] b;

    analitycznie.close();
    numerow.close();

    return blad[naj];
}

double trzypunktowa(double h, int N) {
    double *l, *d, *u, *b, *x, *blad, x0 = A, xn = A;
    fstream konwencjonalnie;
    konwencjonalnie.open("wyniki_trzypunktowe.txt", ios_base::app);
    konwencjonalnie << scientific;
    cout.precision(10);
    l = new double[N];
    d = new double[N];
    u = new double[N];
    b = new double[N];
    x = new double[N];
    blad = new double[N];

    // Z warunku brzegowego
    u[0] = ALFA / h;
    d[0] = BETA - ALFA / h;
    b[0] = -GAMMA;

    // Wyznaczenie srodkowych wyrazow w zapisie wektorowym
    for (int i = 1; i < N - 1; i++) {
        l[i - 1] = P / (h * h) - Q / (2.0 * h);
        d[i] = (-2.0 * P) / (h * h) + R;
        u[i] = P / (h * h) + Q / (2.0 * h);
        b[i] = (x0 + i * h);
    }

    // Z warunku brzegowego
    l[N - 2] = -FI / h;
    d[N - 1] = -FI / h + PSI;
    b[N - 1] = -TETA;

    thomas(l, d, u, b, x, N);

    for (int i = 0; i < N; i++) {
        blad[i] = fabs(x[i] - rozwiazanie_analityczne(xn));
        xn += h;
    }

    int naj = norma_maksimum(blad, N);
    if (N == 1002) {
        for (int i = 0; i < N; i++) {
        konwencjonalnie << x0 << "\t" << x[i] << "\n";
        x0 += h;
        }
    }

    delete[] l;
    delete[] d;
    delete[] u;
    delete[] x;
    delete[] b;

    return blad[naj];
}

int main() {
    double h, blad_num, blad_konw;
    int N; //ilość iteracji

    fstream bledy, rzedy;
    bledy.open("bledy.txt", fstream::out);
    bledy << scientific;
    cout.precision(10);

    for (N = 2; N < 1000000; N += 100) {
        h = (B - A) / (N - 1);
        blad_konw = trzypunktowa(h, N);
        blad_num = numerow(h, N);
        bledy << log10(h) << "\t" << log10(blad_konw) << "\t" << log10(blad_num) << "\n";
    }
    bledy.close();
    return 0;
}
