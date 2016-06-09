#ifndef INTEGRAL_H_INCLUDED
#define INTEGRAL_H_INCLUDED

#include <cmath>
#include <cstdio>
#include <cstdlib>
using namespace std;

class Integral {

public:
    double f(double x); //  Integrated function
    int qpow(int x, int b);

    int stepCal(double a, double b, double res, double eps);
    double trapezoidal(double a, double b, int n);
    double simpson(double a, double b, int n);
    double romberg(double a, double b, double eps);

};

#endif // INTEGRAL_H_INCLUDED
