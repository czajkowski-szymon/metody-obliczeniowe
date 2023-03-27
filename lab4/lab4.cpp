#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

#define MAX_N 50
const double TOLX = 1e-8;
const double TOLF = 1e-8;

double max(double x, double y, double z) {
    double max = x;

    if (y > max) {
        max = y;
    }
    if (z > max) {
        max = z;
    }

    return max;
}

double funkcja_1(double x, double y, double z) {
    return x*x + y*y + z*z - 2.0;
}

double funkcja_2(double x, double y, double z) {
    return x*x + y*y - 1.0;
}

double funkcja_3(double x, double y, double z) {
    return x*x - y;
}

double delta_1(double x, double y, double z) {
    return (x*x - y*y + 2.0*x*x*y - 1.0) / (4.0*x*y + 2.0*x);
}

double delta_2(double x, double y, double z) {
    return (y*y + y - 1.0) / (2.0*y + 1.0);
}

double delta_3(double x, double y, double z) {
    return (z - 1.0) / (2.0*z);
}

void newton(double x, double y, double z) {
    double estymator_bledu = 0, residuum = 0;
    //double wektor_n[3] = {x, y, z};
    double wektor_n1[3] = {0.0, 0.0, 0.0};
    double delta[3] = {0.0, 0.0, 0.0};

    printf("NR ======= X ============ Y ============= Z =========== RESIDUUM =========== BLAD =====\n");

    for (int i = 0; i < MAX_N; i++) {
        delta[0] = delta_1(x, y, z);
        delta[1] = delta_2(x, y, z);
        delta[2] = delta_3(x, y, z);

        wektor_n1[0] = x - delta[0]; 
        wektor_n1[1] = y - delta[1];
        wektor_n1[2] = z - delta[2];

        estymator_bledu = max(fabs(wektor_n1[0] - x), 
                              fabs(wektor_n1[1] - y),
                              fabs(wektor_n1[2] - z));

        residuum = max(fabs(funkcja_1(wektor_n1[0], wektor_n1[1], wektor_n1[2])),
                       fabs(funkcja_2(wektor_n1[0], wektor_n1[1], wektor_n1[2])),
                       fabs(funkcja_3(wektor_n1[0], wektor_n1[1], wektor_n1[2])));

        printf("%2d. %+.10f | %+.10f | %+.10f | %.10e | %.10e\n", i+1, wektor_n1[0], wektor_n1[1], wektor_n1[2], residuum, estymator_bledu);
    
        if (fabs(estymator_bledu) <= TOLX || fabs(residuum) <= TOLF) {
            break;
        }

        x = wektor_n1[0];
        y = wektor_n1[1];
        z = wektor_n1[2];
    }
}

int main() {
    newton(1.0, 2.0, 5.0);

    return 0;
}