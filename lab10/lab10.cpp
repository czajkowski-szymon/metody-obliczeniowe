#include <cmath> 
#include <iostream> 
#include <fstream> 

using namespace std;

double rozwiazanie_analityczne(double t) {
    return 1.0 - exp(-10.0 * (t + atan(t)));
}

double f(double t) {
    return (10.0 * t * t + 20.0) / (t * t + 1.0);
}

double BME(double dt, double tk, double yk) {
    return yk + dt * (f(tk) * (1 - yk));
}


double PME(double dt, double tk, double yk) {
    double val = dt * f(tk + dt);
    return (yk + val) / (1 + val);
}

double PMT(double dt, double tk, double yk) {
    double f_tk = f(tk);
    double f_tk_dt = f(tk + dt);
    return ((-dt / 2.0) * ((yk - 1) * f_tk - f_tk_dt) + yk) / ((dt / 2.0) * f_tk_dt + 1);
}

double bledyBME(double dt) {
    double y = 0.0, yk, max_blad = 0.0, blad, tk = 0.0;
    
    while (tk < 1.0) {
        yk = BME(dt, tk, y);
        blad = fabs(rozwiazanie_analityczne(tk) - y);
        if (blad > max_blad) {
            max_blad = blad;
        }
        tk += dt;
        y = yk;
    }

    return max_blad;
}

double bledyPME(double dt) {
    double y = 0.0, yk, max_blad = 0.0, blad, tk = 0.0;
    
    while (tk < 1.0) {
        yk = PME(dt, tk, y);
        blad = fabs(rozwiazanie_analityczne(tk) - y);
        if (blad > max_blad) {
            max_blad = blad;
        }
        tk += dt;
        y = yk;
    }
    
    return max_blad;
}

double bledyPMT(double dt) {
    double y = 0.0, yk, max_blad = 0.0, blad, tk = 0.0;
    
    while (tk < 1.0) {
        yk = PMT(dt, tk, y);
        blad = fabs(rozwiazanie_analityczne(tk) - y);
        if (blad > max_blad) {
            max_blad = blad;
        }
        tk += dt;
        y = yk;
    }
    
    return max_blad;
}

int main() {
    double blad_bme = 0.0, blad_pme = 0.0, blad_pmt = 0.0;
    ofstream bledy, stabilne, niestabilne;

    bledy.open("bledy.txt");
    stabilne.open("stabilne.txt");
    niestabilne.open("niestabilne.txt");

    double dt = 0.01;
    while (dt > 1e-11) {
        blad_bme = log10(bledyBME(dt));
        blad_pme = log10(bledyPME(dt));
        blad_pmt = log10(bledyPMT(dt));
        bledy << log10(dt) << "\t" << blad_bme << "\t" << blad_pme << "\t" << blad_pmt << endl;
        dt = dt / 2;
    }

    dt = 0.005;
    double wynik_bme = 0.0, wynik_pme = 0.0, wynik_pmt = 0.0;
    double tk = 0.0;
    while (tk < 1.0) {
        wynik_bme = BME(dt, tk, wynik_bme);
        wynik_pme = PME(dt, tk, wynik_pme);
        wynik_pmt = PMT(dt, tk, wynik_pmt);
        stabilne << tk << "\t" << rozwiazanie_analityczne(tk) << "\t" << wynik_bme << "\t" << wynik_pme << "\t" << wynik_pmt << endl;
        tk += dt;
    }
    
    dt = 0.25;
    tk = 0;
    wynik_bme = 0.0;
    while (tk < 5.0) {
        wynik_bme = BME(dt, tk, wynik_bme);
        niestabilne << tk << "\t" << rozwiazanie_analityczne(tk) << "\t" << wynik_bme << endl;
        tk += dt;
    }

    bledy.close();
    stabilne.close();
    niestabilne.close();
    
    return 0;
}