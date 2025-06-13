MN13 – Biblioteka metod numerycznych

Projekt zawiera implementację podstawowych metod numerycznych w języku C++ w postaci biblioteki statycznej. Biblioteka umożliwia rozwiązywanie układów równań liniowych, interpolację, całkowanie numeryczne, rozwiązywanie równań różniczkowych oraz znajdowanie pierwiastków funkcji.

----------------------------Zawartość-------------------------------

naglowek.hpp – interfejs biblioteki

source.cpp – implementacja metod

tests.cpp – testy jednostkowe z użyciem Catch2

main.cpp – przykład użycia w programie głównym (Gauss)

examples - przykłady użycia 

CMakeLists.txt – konfiguracja projektu

catch.hpp – nagłówek frameworka testowego Catch2

------------------Kompilacja i uruchomienie---------------------------

Wymagania: g++, cmake

mkdir build && cd build
cmake ..
make

Uruchomienie programu głównego:

./main

Uruchomienie testów:

./tests
# lub
ctest

-------------------------Przykłady użycia----------------------------

Gauss – rozwiązywanie układu równań:

std::vector<double> A = {2,1,-1, -3,-1,2, -2,1,2};
std::vector<double> b = {8,-11,-3};
auto x = mynum::gauss(A,b,3);

Interpolacja Lagrange’a:

std::vector<double> xi = {0, 1};
std::vector<double> fi = {5, 7};
double y = mynum::lagrange(0.5, xi, fi); // powinno dać 6

Całkowanie – metoda Simpsona:

auto f = [](double x) { return x*x; };
double integral = mynum::simpson(f, 0.0, 1.0, 100);

Bisekcja:

auto f = [](double x) { return std::sin(x); };
int iter;
double root = mynum::bisection(f, 3, 4, 1e-8, 100, iter); // zbliżone do pi

Metoda Heuna (ODE):

auto f = [](double, double y) { return y; };
auto ys = mynum::ode_heun(f, 0.0, 1.0, 0.1, 10);

-----------------------------Licencja--------------------------------

Projekt przeznaczony do celów edukacyjnych w ramach zajęć z metod numerycznych.

--------------------------------------------------------------------
Autor: Paweł Juras
