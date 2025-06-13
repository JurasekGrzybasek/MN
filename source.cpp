// source.cpp
#include "naglowek.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <limits>

#define EPS 1e-10
namespace mynum {

/* ================== Gauss ================== */

/**
 * Rozwiązuje układ równań liniowych Ax = b metodą eliminacji Gaussa.
 * @Ain - macierz współczynników A (w postaci jednowymiarowej tablicy n*n)
 * @bin - wektor wyrazów wolnych b
 * @n - liczba równań (i zmiennych)
 * @return Wektor rozwiązań x (std::vector<double>) w kolejności zmiennych.
 * @throws std::runtime_error jeśli macierz jest osobliwa
 * @throws std::invalid_argument jeśli rozmiary wejścia są nieprawidłowe
 */
std::vector<double> gauss(const std::vector<double>& Ain,
                          const std::vector<double>& bin,
                          int n)
{
    if ((int)Ain.size() != n*n || (int)bin.size() != n)
        throw std::invalid_argument("Invalid system size for Gauss");

    std::vector<double> M = Ain;
    std::vector<double> rhs = bin;

    for (int k=0;k<n;++k) {
        int piv=k; double maxv=std::fabs(M[k*n+k]);
        for (int i=k+1;i<n;++i) {
            double v=std::fabs(M[i*n+k]);
            if (v>maxv){maxv=v;piv=i;}
        }
        if (maxv<EPS) throw std::runtime_error("singular matrix");

        if (piv!=k){
            for(int j=0;j<n;++j) std::swap(M[k*n+j],M[piv*n+j]);
            std::swap(rhs[k],rhs[piv]);
        }

        double diag=M[k*n+k];
        for(int j=k;j<n;++j) M[k*n+j]/=diag;
        rhs[k]/=diag;

        for(int i=0;i<n;++i){
            if(i==k)continue;
            double fct=M[i*n+k];
            for(int j=k;j<n;++j) M[i*n+j]-=fct*M[k*n+j];
            rhs[i]-=fct*rhs[k];
        }
    }
    return rhs;
}

/* ================== Lagrange ================== */

/**
 * Oblicza wartość funkcji interpolacyjnej Lagrange’a w punkcie x.
 * @x - punkt, w którym liczymy wartość
 * @xi - wektor punktów węzłowych
 * @fi - wartości funkcji w punktach xi
 * @return Wartość interpolowaną w punkcie x (double).
 * @throws std::invalid_argument jeśli xi i fi mają różne rozmiary
 */
double lagrange(double x,
                const std::vector<double>& xi,
                const std::vector<double>& fi)
{
    if (xi.size() != fi.size())
        throw std::invalid_argument("xi and fi size mismatch in Lagrange");

    int n=static_cast<int>(xi.size());
    double sum=0.0;
    for(int i=0;i<n;++i){
        double li=1.0;
        for(int j=0;j<n;++j)
            if(j!=i) li*=(x-xi[j])/(xi[i]-xi[j]);
        sum+=fi[i]*li;
    }
    return sum;
}

/* ================== Simpson ================== */

/**
 * Oblicza przybliżoną całkę oznaczoną funkcji f na przedziale [a, b] metodą Simpsona.
 * @f - funkcja podcałkowa
 * @a - początek przedziału
 * @b - koniec przedziału
 * @K - liczba podprzedziałów (musi być parzysta; domyślnie 200)
 * @return Przybliżona wartość całki (double). Zwraca NaN przy niepoprawnych danych.
 */
double simpson(const std::function<double(double)>& f,
               double a,double b,int K)
{
    if (a == b || K < 2)
        return std::numeric_limits<double>::quiet_NaN();

    if(K&1) ++K;
    double h=(b-a)/K;
    double s=f(a)+f(b);
    for(int i=1;i<K;++i)
        s+=(i&1?4.0:2.0)*f(a+i*h);
    return s*h/3.0;
}

/* ================== Heun ================== */

/**
 * Oblicza jeden krok metody Heuna (jawnej) dla równania różniczkowego y' = f(t,y).
 * @f - funkcja pochodnej f(t,y)
 * @t - aktualna wartość t
 * @y - aktualna wartość y
 * @h - krok czasowy
 * @return Przybliżona wartość y po jednym kroku (double)
 */
double step_heun(const std::function<double(double,double)>& f,
                 double t,double y,double h)
{
    double k1=f(t,y);
    double k2=f(t+h,y+h*k1);
    return y+0.5*h*(k1+k2);
}

/**
 * Rozwiązuje równanie różniczkowe y'=f(t,y) metodą Heuna dla n kroków.
 * @f - funkcja pochodnej
 * @t0 - początkowy czas
 * @y0 - początkowa wartość
 * @h - krok czasowy
 * @pn - liczba kroków
 * @return Wektor wartości y w kolejnych krokach (std::vector<double>), łącznie z y0
 */
std::vector<double> ode_heun(const std::function<double(double,double)>& f,
                             double t0,double y0,double h,int n)
{
    std::vector<double> ys; ys.reserve(n+1);
    double t=t0,y=y0; ys.push_back(y);
    for(int i=0;i<n;++i){ y=step_heun(f,t,y,h); t+=h; ys.push_back(y); }
    return ys;
}

/* ================== Bisekcja ================== */

/**
 * Znajduje pierwiastek funkcji f w przedziale [a,b] metodą bisekcji.
 * @f - funkcja, dla której szukamy miejsca zerowego
 * @a - lewy koniec przedziału
 * @b - prawy koniec przedziału
 * @eps - tolerancja błędu
 * @ maxIter - maksymalna liczba iteracji
 * @it - liczba wykonanych iteracji (wyjściowo)
 * @return Wartość przybliżonego pierwiastka (double), lub NaN jeśli warunki są niespełnione.
 */
double bisection(const std::function<double(double)>& f,
                 double a,double b,double eps,int maxIter,int& it)
{
    if(std::isnan(f(a))||std::isnan(f(b))||f(a)*f(b)>0){
        it = 0;
        return std::numeric_limits<double>::quiet_NaN();
    }

    it=0;
    while((b-a)>=eps && it<maxIter){
        double c=0.5*(a+b), fc=f(c);
        if(std::isnan(fc)) return std::numeric_limits<double>::quiet_NaN();
        if(std::fabs(fc)<eps) return c;
        (f(a)*fc<0)? b=c : a=c;
        ++it;
    }
    return 0.5*(a+b);
}

} // namespace mynum
