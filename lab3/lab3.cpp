#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

#define MAX_N 50
const double TOLX = 1e-10;
const double TOLF = 1e-10;

typedef double(*funkcja)(double);

// sin^2(x/4) - x = 0

double funkcja_1(double x) {
    return sin(x / 4.0) * sin(x / 4.0) - x;
}

double funkcja_pochodna_1(double x) {
    return 0.25 * sin(x / 2.0) - 1.0;
}

double fi_1(double x) {
    return sin(x / 4.0) * sin(x / 4.0);
}

// tan(2x) - x - 1 = 0

double funkcja_2(double x) {
    return tan(2.0 * x) - x - 1.0;
}

double fi_2(double x) {
    return 0.5 * atan(x + 1.0);
}

double funkcja_pochodna_2(double x) {
    return -1.0 + 2.0 / (cos(2.0 * x) * cos(2.0 * x));
}

void picard(funkcja funkcja_poczatkowa, funkcja fi, double x) {
    printf("PICARD ------------------------------------------------\n");
    printf("NR ====== X =========== RESIDUUM =========== BLAD =====\n");
    
    double estymator_bledu, residuum, x_n;

    for (int i = 0; i < MAX_N; i++) {
        // obliczenie przybliżenia
        x_n = fi(x);
        
        // obliczenie estymatora bledu
        estymator_bledu = fabs(x_n - x);
        
        // obliczenie residuum
        residuum = fabs(funkcja_poczatkowa(x_n));
        
        printf("%2d. %+.10f | %.10e | %.10e\n", i+1, x_n, residuum, estymator_bledu);

        // sprawdzenie warunków stopu
        if (estymator_bledu <= TOLX || residuum <= TOLF) {
            break;
        }

        x = x_n;
    }
}

void bisekcja(funkcja funkcja_poczatkowa, double a, double b) {
    printf("BISEKCJA ----------------------------------------------\n");
    printf("NR ====== X =========== RESIDUUM =========== BLAD =====\n");
    
    double estymator_bledu, residuum, x;

    for (int i = 0; i < MAX_N; i++) {
        // obliczenie przybizenia
        x = (a + b) / 2.0;
        
        // obliczenie estymatora bledu
        estymator_bledu = fabs((b - a) / 2.0);
        
        // residuum
        residuum = fabs(funkcja_poczatkowa(x));

        printf("%2d. %+.10f | %.10e | %.10e\n", i+1, x, residuum, estymator_bledu);

        // wybor nowego przedzialu
        if (funkcja_poczatkowa(a) < 0 && funkcja_poczatkowa(x) > 0 || funkcja_poczatkowa(a) > 0 && funkcja_poczatkowa(x) < 0) {
            b = x;
        } else {
            a = x;
        }

        // sprawdzenie warunków stopu
        if (estymator_bledu <= TOLX || residuum <= TOLF) {
            break;
        }
    }
}

void newton(funkcja funkcja_poczatkowa, funkcja funkcja_pochodna, double x) {
    printf("NEWTON ------------------------------------------------\n");
    printf("NR ====== X =========== RESIDUUM =========== BLAD =====\n");
    
    double estymator_bledu, residuum, x_n; 
    
    for (int i = 0; i < MAX_N; i++) {
        // obliczenie przyblizenia
        x_n = x - (funkcja_poczatkowa(x) / funkcja_pochodna(x));
        
        // obliczenie estymatora bledu
        estymator_bledu = fabs(x_n - x);
        

        // obliczenie residuum
        residuum = fabs(funkcja_poczatkowa(x_n));

        printf("%2d. %+.10f | %.10e | %.10e\n", i+1, x_n, residuum, estymator_bledu);
    
        // sprawdzenie warunków stopu
        if (estymator_bledu <= TOLX || residuum <= TOLF) {
            break;
        }

        x = x_n; 
    }
}

void sieczne(funkcja funkcja_poczatkowa, double x1, double x2) {
    printf("SIECZNE -----------------------------------------------\n");
    printf("NR ====== X =========== RESIDUUM =========== BLAD =====\n");

    double estymator_bledu, residuum, x3; 

    for (int i = 0; i < MAX_N; i++) {
        // obliczenie porzyblizenia
        x3 = x2 - funkcja_poczatkowa(x2) / ((funkcja_poczatkowa(x2) - funkcja_poczatkowa(x1)) / (x2 - x1));

        // obliczenie estymatora bledu
        estymator_bledu = fabs(x3 - x2);

        // obliczenie residuum
        residuum = fabs(funkcja_poczatkowa(x3));
    
        printf("%2d. %+.10f | %.10e | %.10e\n", i+1, x3, residuum, estymator_bledu);

        // sprawdzenie warunku stopu
        if (estymator_bledu <= TOLX || residuum <= TOLF) {
            break;
        }

        x1 = x2;
        x2 = x3;
    }
}

int main() {
    cout << "sin^2(x/4) - x = 0" << endl;

    picard(funkcja_1, fi_1, 1.0);
    bisekcja(funkcja_1, -0.7, 1.2);
    newton(funkcja_1, funkcja_pochodna_1, 1.0);
    sieczne(funkcja_1, 0.4, 1.2);

    cout << endl;
    cout << "tan(2x) - x - 1 = 0" << endl;

    picard(funkcja_2, fi_2, 1.0);
    bisekcja(funkcja_2, -0.7, 1.2);
    newton(funkcja_2, funkcja_pochodna_2, 0.4);
    sieczne(funkcja_2, 0.7, 1.3);

    return 0;
}