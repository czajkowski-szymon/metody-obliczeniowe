#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

using namespace std;

#define N 100

double funkcja(double x) {
    return (1.0 - exp(-x)) / x;
}

double szereg_taylora(double x) {
    double wynik = 1.0;
    double znak = -1.0;
    double krok = 1.0;

    for (int i = 1; i < N; i++){
        krok *= (-1.0) * (x / (i + 1));
        wynik += krok;
        znak = -znak;
    }

    return wynik;
}

int main() {
    ifstream plik_dane;
    ofstream plik_wynik;
    ofstream plik_wynik2;
    plik_dane.open("dane.txt");
    plik_wynik.open("wynik.txt");
    plik_wynik2.open("wynik2.txt");
    double log_10, x, wartosc_z_pliku, wartosc_obliczona, blad, blad_log;

    while (!plik_dane.eof()) {
        plik_dane >> log_10;
        plik_dane >> x;
        plik_dane >> wartosc_z_pliku;

        wartosc_obliczona = funkcja(x);
        blad = abs((wartosc_obliczona - wartosc_z_pliku) / wartosc_z_pliku);
        blad_log = log10(blad);
        
        plik_wynik << setprecision(20) << scientific << log_10 << " " << blad_log << endl;

        if (log_10 < 0.0) {
            wartosc_obliczona = szereg_taylora(x);
            blad = abs((wartosc_obliczona - wartosc_z_pliku) / wartosc_z_pliku); 
            blad_log = log10(blad);
            plik_wynik2 << setprecision(20) << scientific << log_10 << " " << blad_log << endl; 
        } else {
            wartosc_obliczona = funkcja(x);
            blad = abs((wartosc_obliczona - wartosc_z_pliku) / wartosc_z_pliku); 
            blad_log = log10(blad);
            plik_wynik2 << setprecision(20) << scientific << log_10 << " " << blad_log << endl;
        }
    }
    
    plik_dane.close();
    plik_wynik.close();

    return 0;
}