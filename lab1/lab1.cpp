#include <iostream>

int main() {
    int bitDouble = 0;
    double precDouble = 1.0;
    double tmpDouble = 1.1;

    while (tmpDouble > 1.0) {
        precDouble /= 2.0;
        tmpDouble = 1.0 + precDouble;
        if (tmpDouble != 1.0) {
            bitDouble++;
        }
    }

    double epsilonDouble = 2.0 * precDouble;
    std::cout << "double: " <<  bitDouble << " " << epsilonDouble << std::endl;

    int bitFloat = 0;
    float precFloat = 1.0f;
    float tmpFloat = 1.1f;

    while (tmpFloat > 1.0f) {
        precFloat /= 2.0f;
        tmpFloat = 1.0f + precFloat;
        if (tmpFloat != 1.0f) {
            bitFloat++;
        }
    }

    float epsilonFloat = 2.0f * precFloat;
    std::cout << "float: " << bitFloat << " " << epsilonFloat << std::endl;

    return 0;
}

