#include "Integral.h"

const int MAXN = 10000;

double Integral::f(double x) {
    return sqrt(4 - sin(x) * sin(x));
}

int Integral::qpow(int x, int b) {
    int res = 1;
    while(b) {
        if(b & 1) res *= x;
        x *= x;
        b /= 2;
    }
    return res;
}

int Integral::stepCal(double a, double b, double res, double eps){
    int l = 1, r = 100;
    while(l <= r){
        int mid = (l + r) >> 1;
        if(abs(trapezoidal(a, b, mid) - res) >= eps) l = mid + 1;
        else r = mid - 1;
    }
    printf("%.10f\nstep == %d\n",trapezoidal(a, b, l), l);
    return l;
}

double Integral::trapezoidal(double a, double b, int n){
    double sum = 0, h = (b - a) / n;
    sum += f(a) + f(b);
    for(int i = 1; i < n ; i ++) sum += 2 * f(a + i * h);
    return h * sum / 2;
}

double Integral::simpson(double a, double b, int n) {

    //simpson: f(x) integral(a -> b)  = (b - a) * (f(a) + 4 * f((a+b)/2) + f(b) / 6

    double sum1, sum2, sum3, x, h;
    sum1 = sum2 = sum3 = 0;

    h = (b - a)/(2*n);
    sum1 = f(a) + f(b);
    for(int i = 1; i < 2 * n; i ++) {
        x = a + i * h;
        if(i & 1) sum2 += f(x);
        else sum3 += f(x);
    }
    return h * (sum1 + 4 * sum2 + 2*sum3) / 3;
}

double Integral::romberg(double a, double b, double eps) {
    double T[2][MAXN];
    int k = 0, now = 0, pre = 1 - now;
    T[now][0] = (b - a) * (f(a) + f(b)) / 2;

    do {
        k ++;
        now = 1 - now;
        pre = 1 - now;

        int NUM = qpow(2, k - 1);
        T[now][0] = 0;
        for(int i = 1; i <= NUM; i ++)
            T[now][0] += f(a + (2 * i - 1) * (b - a) / qpow(2, k));
        T[now][0] *= (b - a) / qpow(2, k);
        T[now][0] += T[pre][0] / 2;

        for(int i = 1; i <= k; i ++)
            T[now][i] = (qpow(4, i) * T[now][i-1] - T[pre][i-1]) / (qpow(4, i) - 1);

    } while(abs(T[now][k - 1] - T[pre][k - 1] >= eps));
    return T[now][k - 1];
}
