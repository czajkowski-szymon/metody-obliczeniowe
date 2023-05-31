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

const double H = 0.1;

const double LAMBDA = 1.0;

void obroc_macierz(double **macierz, int N, int M)
{
    for (int i = 0; i < N / 2; i++)
    {
        for (int j = 0; j < M; j++)
        {
            swap(macierz[i][j], macierz[N - 1 - i][j]);
        }
    }
}

double **stworz_macierz(int N, int M)
{
    double **macierz = new double *[N];
    for (int i = 0; i < N; i++)
    {
        macierz[i] = new double[M];
    }
    return macierz;
}

void usun_macierz(double **macierz, int N)
{
    for (int i = 0; i < N; i++)
    {
        delete[] macierz[i];
    }
    delete[] macierz;
}

double *siatka_kroku(double h, int M)
{
    double *wektor = new double[M];
    double x = POCZATEK_PRZEDZIALU;

    for (int i = 0; i < M; i++)
    {
        wektor[i] = x;
        x += h;
    }

    return wektor;
}

void zapisz_do_pliku_w(double *wektor, int N, string nazwa)
{
    ofstream plik;
    plik.open(nazwa);
    double *kroki = siatka_kroku(H, N);
    for (int i = 0; i < N; i++)
    {
        plik << kroki[i] << "\t" << wektor[i] << endl;
    }
}

void zapisz_do_pliku_m(double **macierz, int N, int M, string nazwa)
{
    ofstream plik;
    plik.open(nazwa);
    double *kroki = siatka_kroku(H, N);
    for (int i = 0; i < N; i++)
    {
        plik << kroki[i] << "\t";
        for (int j = 0; j < M; j++)
        {
            plik << macierz[i][j] << "\t";
        }
        plik << endl;
    }

    plik.close();
}

double rozwiazanie_analityczne(double x, double t)
{
    return 0.5 * calerfpack::erfc_l(x / 2 * sqrt(D * t));
}

double **oblicz_rozwiazanie_analityczne(double h, double dt, int N, int M)
{
    double **macierz_rozwiazan = stworz_macierz(N, M);
    double x = POCZATEK_PRZEDZIALU;
    double t = T_MIN;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            macierz_rozwiazan[i][j] = rozwiazanie_analityczne(x, t);
            x += h;
        }
        x = POCZATEK_PRZEDZIALU;
        t += dt;
    }

    return macierz_rozwiazan;
}

void warunek_poczatkowy(double **macierz, double h, int M)
{
    double x = POCZATEK_PRZEDZIALU;
    for (int i = 0; i < M; i++)
    {
        if (x < 0.0)
        {
            macierz[0][i] = 1.0;
        }
        else
        {
            macierz[0][i] = 0.0;
        }
        x += h;
    }
}

void warunek_brzegowy(double **macierz, int N, int M)
{
    for (int i = 0; i < N; i++)
    {
        macierz[i][0] = 1.0;
        macierz[i][M - 1] = 0.0;
    }
}

double norma_maksimum(double *wektor, int N)
{
    double max = fabs(wektor[0]);
    for (int i = 1; i < N; i++)
    {
        if (max < fabs(wektor[i]))
        {
            max = wektor[i];
        }
    }
    return max;
}

double* policz_bledy(double **analityczne, double **numeryczne, int n, int m)
{
    double **bledy = stworz_macierz(n, m);
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            bledy[i][j] = fabs(numeryczne[i][j] - analityczne[i][j]);
        }
    }
    double* max_blad = new double[n];
    for (int i = 0; i < n; i++)
    {
        max_blad[i] = norma_maksimum(bledy[i], m);
    }
    usun_macierz(bledy, n);
    return max_blad;
}

double **transponuj(double **macierz, int N, int M)
{
    double **result = stworz_macierz(M, N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < M; ++j)
        {
            result[j][i] = macierz[i][j];
        }
    }

    return result;
}

void thomas(double **a, double *b, double *x, int m)
{
    double *l = new double[m];
    double *d = new double[m];
    double *u = new double[m];

    d[0] = a[0][0];

    for (int i = 1; i < m - 1; i++)
    {
        u[i] = a[i][i + 1];
        d[i] = a[i][i];
        l[i - 1] = a[i][i - 1];
    }

    d[m - 1] = a[m - 1][m - 1];

    for (int i = 1; i < m; i++)
    {
        d[i] = d[i] - ((l[i - 1] / d[i - 1]) * u[i - 1]);
        b[i] = b[i] - ((l[i - 1] / d[i - 1]) * b[i - 1]);
    }

    x[m - 1] = b[m - 1] / d[m - 1];

    for (int i = m - 2; i >= 0; i--)
    {
        x[i] = (b[i] - u[i] * x[i + 1]) / d[i];
    }
}

double **laasonen_thomas_pom(double h, int n, int m)
{
    double **rozwiazanie = stworz_macierz(n, m);
    warunek_poczatkowy(rozwiazanie, h, m);
    warunek_brzegowy(rozwiazanie, n, m);
    double *b = new double[m];
    double *wynik = new double[m];

    for (int i = 0; i < m; i++)
    {
        b[i] = 0.0;
        wynik[i] = 0.0;
    }

    double **a = stworz_macierz(m, m);

    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < m; j++)
        {
            a[i][j] = 0.0;
        }
    }

    for (int k = 1; k < n; k++)
    {
        a[0][0] = 1.0;
        b[0] = rozwiazanie[k - 1][0];

        for (int i = 1; i < m - 1; i++)
        {
            a[i][i + 1] = LAMBDA;
            a[i][i] = -(1.0 + 2.0 * LAMBDA);
            a[i][i - 1] = LAMBDA;
            b[i] = -rozwiazanie[k - 1][i];
        }

        b[m - 1] = 0.0;
        a[m - 1][m - 1] = 1.0;

        thomas(a, b, wynik, m);

        for (int i = 1; i < m - 1; i++)
        {
            rozwiazanie[k][i] = wynik[i];
        }
    }

    usun_macierz(a, n);
    delete[] b;
    delete[] wynik;
    return rozwiazanie;
}

void laasonen_thomas()
{
    double h = 0.1;
    double dt = (LAMBDA * h * h) / D;
    int n = (T_MAX - T_MIN) / dt;
    int m = (KONIEC_PRZEDZIALU - POCZATEK_PRZEDZIALU) / h;
    double** analityczne = oblicz_rozwiazanie_analityczne(h, dt, n, m);
    double** numeryczne = laasonen_thomas_pom(h, n, m);
    obroc_macierz(numeryczne, n, m);
    zapisz_do_pliku_m(transponuj(analityczne, n, m), m, n, "rozwiazanie_analityczne.txt");
    zapisz_do_pliku_m(transponuj(numeryczne, n, m), m, n, "lt_numerycznie.txt");
    usun_macierz(analityczne, n);
    usun_macierz(numeryczne, n);
}

void laasonen_thomas_bledy()
{
    double kroki[4] = {0.1, 0.01, 0.005, 0.001};
    double h = 0.1;
    double dt;
    int n, m;
    ofstream plik;
    //plik.open("lt_bledy.txt");
    for (h; h > 1e-4; h /= 2) {
        // double h = kroki[4];
        // cout << h << endl;
        dt = (LAMBDA * h * h) / D;
        n = (T_MAX - T_MIN) / dt;
        m = (KONIEC_PRZEDZIALU - POCZATEK_PRZEDZIALU) / h;
        double** analityczne = oblicz_rozwiazanie_analityczne(h, dt, n, m);
        double** numeryczne = laasonen_thomas_pom(h, n, m);
        double* bledy = policz_bledy(analityczne, numeryczne, n, m);
        usun_macierz(numeryczne, n);
        usun_macierz(analityczne, n);
        //plik << log10(h) << "\t" << log10(bledy[n - 1]) << endl;
        cout << log10(h) << "\t" << log10(bledy[n - 1]) << endl;
        delete[] bledy;
    }

    plik.close();
}

int main()
{
    /*
     *  LAASONEN - THOMAS
     */

    //laasonen_thomas();
    laasonen_thomas_bledy();


    // ofstream bledy_plik;
    // bledy_plik.open("lt_bledy.txt");
    // for (h; h > 1e-2; h /= 2) {
    //     double dt = (LAMBDA*h*h) / D;
    //     int N = (T_MAX - T_MIN) / dt;
    //     int M = (KONIEC_PRZEDZIALU - POCZATEK_PRZEDZIALU) / h;
    //     rozwiazanie_analityczne = oblicz_rozwiazanie_analityczne(h, dt, N, M);
    //     laasonen_thomas_sol = laasonen_thomas(h, N, M);
    //     bledy = policz_bledy(rozwiazanie_analityczne, laasonen_thomas_sol, N, M);
    //     maks_blad = norma_maksimum(bledy[N-1], M);
    //     bledy_plik << log10(h) << "\t" << log10(maks_blad) << endl;
    // }

    return 0;
}