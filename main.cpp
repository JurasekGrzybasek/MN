#include "naglowek.hpp"
#include <iostream>
#include <cstddef>
#include <vector>


int main() {
    std::vector<double> A = { 2,1,-1,  -3,-1,2,  -2,1,2 };
    std::vector<double> b = { 8,-11,-3 };

    auto x = mynum::gauss(A,b,3);
    std::cout << "x = [" << x[0] << ", " << x[1] << ", " << x[2] << "]\n";
}

