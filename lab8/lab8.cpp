#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

#define N_MAX 30

string nazwa(int i) {
    switch (i)
    {
    case 0:
        return "roznica progresywna dwupunktowa (lewy brzeg) (teoria: 1)";
        break;
    case 1:
        return "roznica centralna dwupunktowa (srodek) (teoria: 2)";
        break;
    case 2:
        return "roznica wsteczna dwupunktowa (prawy brzeg) (teoria: 1)";
        break;
    case 3:
        return "roznica progresywna trzypunktowa (lewy brzeg) (teoria: 2)";
        break;
    case 4:
        return "roznica progresywna dwupunktowa (srodek) (teoria: 1)";
        break;
    case 5:
        return "roznica wsteczna dwupunktowa (srodek) (teoria: 1)";
        break;
    case 6:
        return "roznica wsteczna trzypunktowa (prawy brzeg) (teoria: 2)";
        break;
    default:
        return "-";
        break;
    }
}

template < typename T >
T f(T x) {
    return sin(x);
}

template < typename T >
T pochodna_f(T x) {
    return cos(x);
}

template < typename T >
T roznica_progresywna_2pkt(T x, T h) {
    return (sin(x + h) - sin(x)) / h;
}

template < typename T >
T roznica_wsteczna_2pkt(T x, T h) {
    return (sin(x) - sin(x - h)) / h;
}

template < typename T >
T roznica_centralna_2pkt(T x, T h) {
    return (sin(x + h) - sin(x - h)) / (2.0 * h);
}

template < typename T >
T roznica_progresywna_3pkt(T x, T h) {
    return ((-3.0 / 2.0) * sin(x) + 2.0 * sin(x + h) - (1.0 / 2.0) * sin(x + 2.0 * h)) / h;
}

template < typename T >
T roznica_wsteczna_3pkt(T x, T h) {
    return ((1.0 / 2.0) * sin(x - 2.0*h) - 2.0 * sin(x - h) + (3.0 / 2.0) * sin(x)) / h;
}

template<typename T>
T **oblicz(T** tablica) {
    T** talibca_bledow = new T*[N_MAX];    // Tablica przechowujaca bledy w kolejnych iteracjach
    T h = 0.1;    // Pierwszy krok

    // Przedzial [a, b] gdzie mid to srodek
    T a = 0.0;
    T b = M_PI / 2.0;
    T mid = (a + b) / 2.0;

    for (int i = 0; i < N_MAX; ++i) {
        tablica[i] = new T[8];           // Stworzenie wiersza na wyniki
        talibca_bledow[i] = new T[8];    // Stworzenie wiersza na bledy

        tablica[i][0] = roznica_progresywna_2pkt(a, h);
        talibca_bledow[i][0] = fabs(pochodna_f(a) - tablica[i][0]);

        tablica[i][1] = roznica_centralna_2pkt(mid, h);
        talibca_bledow[i][1] = fabs(pochodna_f(mid) - tablica[i][1]);

        tablica[i][2] = roznica_wsteczna_2pkt(b, h);
        talibca_bledow[i][2] = fabs(pochodna_f(b) - tablica[i][2]);

        tablica[i][3] = roznica_progresywna_3pkt(a, h);
        talibca_bledow[i][3] = fabs(pochodna_f(a) - tablica[i][3]);

        tablica[i][4] = roznica_progresywna_2pkt(mid, h);
        talibca_bledow[i][4] = fabs(pochodna_f(mid) - tablica[i][4]);

        tablica[i][5] = roznica_wsteczna_2pkt(mid, h);
        talibca_bledow[i][5] = fabs(pochodna_f(mid) - tablica[i][5]);

        tablica[i][6] = roznica_wsteczna_3pkt(b, h);
        talibca_bledow[i][6] = fabs(pochodna_f(b) - tablica[i][6]);

        talibca_bledow[i][7] = h;

        h *= 0.1;    
    }

    cout << endl << "Rzedy dokladnosci: " << endl;
    for (int i = 0; i < 7; i++) {
        cout << nazwa(i) << ": ";
        cout << (log10(talibca_bledow[1][i]) - log10(talibca_bledow[0][i])) /
                (log10(talibca_bledow[1][7]) - log10(talibca_bledow[0][7])) <<  endl;
    }
    cout <<  endl;
    return talibca_bledow;
}

template<typename T>
void zapisDoPliku(T** tablica, string nazwa) {
    fstream plik;
    plik << scientific;
    plik.precision(16);
    plik.open(nazwa, ios::out);

    if (plik.good()) {
        for (int i = 0; i < N_MAX; i++) {
        plik << log10(tablica[i][7]) << " ";
        for (int j = 0; j < 7; j++) {
            if (fabs(log10(tablica[i][j])) > 0)
            plik << log10(tablica[i][j]) << " ";
            else
            plik << "0 ";
        }
        plik << "\n";
        }
        plik.close(); 
    } else {
        cout << "Blad otwarcia pliku " << endl;
    }
}

int main() {
    cout << scientific;
    cout.precision(10);

    cout << "FLOAT\n";
    float** pochodna_float = new float*[N_MAX]; 
    pochodna_float = oblicz(pochodna_float);
    zapisDoPliku(pochodna_float, "pochodna_float.txt");

    cout << "DOUBLE\n";
    double** pochodna_double = new double*[N_MAX];
    pochodna_double = oblicz(pochodna_double);
    zapisDoPliku(pochodna_double, "pochodna_double.txt"); 
}