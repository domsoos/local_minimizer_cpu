#include "../src/optimization.h"
#include <cassert>
#include <cmath>
#include <vector>
#include <functional>

// Quadratic function for testing: f(x) = (x-3)^2
double testFunction(std::vector<double> &x) {
    double sum = 0.0;
    for (auto &val : x) {
        sum += std::pow(val - 3.0, 2);
    }
    return sum;
}

std::vector<double> testFunctionDerivative(const std::vector<double> &x) {
    std::vector<double> gradient(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        gradient[i] = 2.0 * (x[i] - 3.0);
    }
    return gradient;
}

void testLineSearch() {
    std::vector<double> x = {0.0}; // Initial point
    std::vector<double> p = {1.0}; // Search direction (towards the minimum)
    std::vector<double> g = testFunctionDerivative(x); // Derivative at initial point

    double f0 = testFunction(x); // Initial function value
    std::cout << "starting point = " << f0 <<std::endl;
    double d0 = directional_derivative(g, p);
    //                                    f      ,  f0         ,    d0,
    double alpha = dlib_line_search(testFunction, testFunction(x),d0,
            0.15,                  // rho
            0.9,                  // sigma
            1e-12,                // min_f
            100,x,p,15);

    // Apply the step to x
    std::vector<double> x_new = { x[0] + alpha * p[0] };

    // Test that the step size is positive and not unreasonably large
    assert(alpha >= 1e-11 && alpha < 10.0);
    // Test that the new function value is less than the original
    std::cout << "predicted function value after step size == " << testFunction(x_new)<<std::endl;
    assert(testFunction(x_new) < f0);
}

int main() {
    testLineSearch();

    return 0;
}
