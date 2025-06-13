#pragma once
#include <cstddef>   
#include <vector>
#include <string>
#include <functional>

namespace mynum {

// ---------- Gauss ----------
bool read_system(const std::string& filename,
                 std::vector<double>& A,
                 std::vector<double>& b,
                 int& n);

std::vector<double> gauss(const std::vector<double>& A,
                          const std::vector<double>& b,
                          int n);

// ---------- Lagrange ----------
double lagrange(double x,
                const std::vector<double>& xi,
                const std::vector<double>& fi);

// ---------- Simpson + LSQ ----------
double simpson(const std::function<double(double)>& f,
               double a, double b, int K = 200);

std::vector<double> poly_fit_lsq(const std::function<double(double)>& f,
                                 int m,
                                 double a, double b,
                                 int K = 200);

// ---------- ODE – Heun ----------
double step_heun(const std::function<double(double,double)>& f,
                 double t, double y, double h);

std::vector<double> ode_heun(const std::function<double(double,double)>& f,
                             double t0, double y0,
                             double h, int n);

// ---------- Bisekcja ----------
double bisection(const std::function<double(double)>& f,
                 double a, double b,
                 double eps, int maxIter,
                 int& iters);

} // namespace mynum

