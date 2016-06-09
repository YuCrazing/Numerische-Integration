#include <cstdio>
#include "Integral.h"
using namespace std;

const double pi = acos(-1.0);

int main(){

    int step = 10;
    double a = 0;
    double b = pi / 4;
    double eps = 0.0000001;
    Integral integral;
    printf("step = %d:\n", step);
    printf("Trapezoidal: %.10f\n", integral.trapezoidal(a, b, step));
    printf("Simpson: %.10f\n", integral.simpson(a, b, step));
    printf("Romberg: %.10f\n", integral.romberg(a, b, eps));
    integral.stepCal(a, b, 1.5343916, 0.00001);

    return 0;
}
