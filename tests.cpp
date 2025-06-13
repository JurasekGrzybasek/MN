#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "naglowek.hpp"
#include <cmath>
#include <limits>

using namespace mynum;

/* --- Gauss --- */
TEST_CASE("Gauss 3x3") {
    std::vector<double> A = {2,1,-1, -3,-1,2, -2,1,2};
    std::vector<double> b = {8,-11,-3};
    auto x = gauss(A,b,3);
    REQUIRE( Approx(2)  == x[0] );
    REQUIRE( Approx(3)  == x[1] );
    REQUIRE( Approx(-1) == x[2] );
}

TEST_CASE("Gauss singular") {
    std::vector<double> A={1,2, 2,4};
    std::vector<double> b={3,6};
    REQUIRE_THROWS_AS( gauss(A,b,2), std::runtime_error );
}

TEST_CASE("Gauss size mismatch") {
    std::vector<double> A={1,2,3,4};
    std::vector<double> b={5};
    REQUIRE_THROWS_AS( gauss(A,b,2), std::invalid_argument );
}

/* --- Lagrange --- */
TEST_CASE("Lagrange linear") {
    std::vector<double> xi={0,1};
    std::vector<double> fi={5,7};
    REQUIRE( lagrange(0.5,xi,fi)==Approx(6) );
}

TEST_CASE("Lagrange bad input") {
    std::vector<double> xi={0,1,2};
    std::vector<double> fi={5,6};
    REQUIRE_THROWS_AS( lagrange(1.0,xi,fi), std::invalid_argument );
}

/* --- Simpson --- */
double f_quad(double x){ return 3*x*x+2*x+1; }

TEST_CASE("Simpson exact quad") {
    REQUIRE( simpson(f_quad,0,1,120)==Approx(3).epsilon(1e-8) );
}

TEST_CASE("Simpson invalid range") {
    REQUIRE( std::isnan(simpson(f_quad, 1.0, 1.0, 20)) );
}

/* --- Heun --- */
TEST_CASE("Heun y'=y") {
    auto f=[](double, double y){return y;};
    auto ys=ode_heun(f,0,1,0.05,20);
    REQUIRE( ys.back()==Approx(std::exp(1)).epsilon(2e-2) );
}

/* --- Bisection --- */
TEST_CASE("Bisection sin") {
    auto f=[](double x){return std::sin(x);};
    int it;
    double root=bisection(f,3,4,1e-10,100,it);
    REQUIRE( root==Approx(M_PI).epsilon(1e-9) );
}

TEST_CASE("Bisection invalid interval") {
    auto f=[](double x){return x*x + 1;}; // no real roots
    int it;
    double r = bisection(f, -1, 1, 1e-6, 100, it);
    REQUIRE( std::isnan(r) );
}

