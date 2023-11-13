#include "../src/test_functions.h"
#include "../lib/utility.h"
#include <cassert>
#include <vector>
#include <cmath>

const double epsilon = 1e-6;

void testPrintVector(const std::vector<double> x0,const double tru) {
    for(auto & val: x0) {
        std::cout << val << " ";
        assert(std::abs(val) == tru);
    }
    std::cout << std::endl;
}

// Function to test the Rosenbrock function
void testRosenbrockFunction() {
    std::vector<double> x0 = {1.0, 1.0};
    std::vector<double> x1 = {-1.2, 1.0};

    double f0 = rosenbrock(x0);
    double f1 = rosenbrock(x1);
    std::cout << "rosenbrock funinit = "<<std::abs(f0) <<std::endl;
    std::cout << "rosenbrock funmin = "<<std::abs(f1) <<std::endl;
    
    assert(std::abs(f1 - 24.2) < epsilon); // absolute margin comparison -- floating point representation
    assert(std::abs(f0) == 0.0); // global min
}

// Function to test the gradient of the Rosenbrock function
void testRosenbrockGradient() {
    std::vector<double> x = {1.0, 1.0}; // Gradient should be zero at the minimum
    std::vector<double> grad = rosenbrock_gradient(x);
    for (auto& val : grad) {
        assert(std::abs(val) == 0); // Each component of the gradient should be close to 0
    }
}

/* Woods Function and its Gradient Tests */
void testWoodsFunction() {
    std::vector<double> x0 = {-3.0, -1.0, -3.0, -1.0};
    std::vector<double> x1 = {1.0,1.0,1.0,1.0};

    double f0 = woods(x0);
    double f1 = woods(x1);

    std::cout << "woods funinit = " << f0 << std::endl;
    std::cout << "woods funmin  = " << f1 << std::endl;

    assert(std::abs(f0 - 19192.0) < epsilon); // init 
    assert(std::abs(f1 - 0.0) < epsilon); // global min
}

void testWoodsGradient() {
    std::vector<double> x0 = {1.0, 1.0, 1.0, 1.0};
    std::vector<double> x1 = {2.0, 2.0, 2.0, 2.0};

    std::vector<double> d0 = woods_derivative(x0);
    std::vector<double> d1 = woods_derivative(x1);

    std::cout << "woods der funmin = ";
    testPrintVector(d0, 0.0);

    assert(std::abs(d1[0] - 1602.0) < epsilon);
    assert(std::abs(d1[1] - -360.0) < epsilon);
    assert(std::abs(d1[2] - 1442.0) < epsilon);
    assert(std::abs(d1[3] - -320.0) < epsilon);
}

/* Powell Quartic and its Gradient Tests*/
void testPowellQuartic() {
    std::vector<double> x0 = {3.0, -1.0, 0.0, 1.0};
    double f0 = powell_quartic(x0);
    std::cout << "powell fun = " << f0 << std::endl;
    assert(std::abs(f0 - 215.0) < epsilon);
}

void testPowellQuarticGradient() {
    std::vector<double> x0 = {0.0, 0.0, 0.0, 0.0}; // funmin = 0
    std::vector<double> x1 = {1.0, 1.0, 1.0, 1.0}; // random point

    std::vector<double> d0 = powell_quartic_derivative(x0);
    std::vector<double> d1 = powell_quartic_derivative(x1);

    std::cout << "powell der funmin = ";
    testPrintVector(d0, 0.0);
    std::cout << "powell der random = ";
    assert(std::abs(d1[0] - 22.0) < epsilon);
    assert(std::abs(d1[1] - 216.0) < epsilon);
    assert(std::abs(d1[2] - 8.0) < epsilon);
    assert(std::abs(d1[3] - 0.0) < epsilon);
}

/* Fletcher Powell Helical Valley Function and its Gradient Tests */
void testFletcherPowellHelicalValley() {
    std::vector<double> x0 = {-1.0, 0.0, 0.0};
    std::vector<double> x1 = {1.0, 0.0, 0.0};

    double f0 = helical_valley(x0);
    double f1 = helical_valley(x1);

    std::cout << "helivalley init = " << f0 << std::endl;
    std::cout << "helivalley min  = " << f1 << std::endl;
    assert(std::abs(f0) == 2500.0);
    assert(std::abs(f1) == 0.0); 
}

void testFletcherPowellHelicalValleyDerivative() {
    std::vector<double> x0 = {1.0, 0.0, 0.0};
    std::vector<double> x1 = {1.0, 1.0, 0.0};

    std::vector<double> d0 = helical_valley_derivative(x0);
    std::vector<double> d1 = helical_valley_derivative(x1);
    std::cout << "helivalley der funmin: ";
    testPrintVector(d0, 0.0);
    std::cout << "helivalley der random: ";
    print_vector(d1);
    assert(std::abs(d1[0] - 0.585786) < epsilon);
    assert(std::abs(d1[1] - 0.585786) < epsilon);
    assert(std::abs(d1[2] - 0.0) < epsilon);
}

/* Fletcher Powell Trigonometric and its Gradient Tests*/
void testFletcherPowellTrigonometric() {
    double x = distr(eng);
    std::vector<double> x0 = {x, x, x};
    double f0 = fletcher_powell_trig(x0);
    std::cout << "fp trig = " << f0 << std::endl;
    assert(std::abs(f0) == 0.0);
}
void testFletcherPowellTrigonometricDerivative() {
    double x = distr(eng);
    std::vector<double> x0 = {x, x, x};
    std::vector<double> d0 = fletcher_powell_trig_derivative(x0);
    std::cout << "fp trig der = ";
    testPrintVector(d0, 0.0);
}

void testThermisterProblem() {
    std::vector<double> x0 = {0.02, 4000, 250};
    double f0 = thermister(x0);
    std::cout << "thermister function = " << f0 << std::endl;
    assert(std::abs(f0 - 1.69e+9) < epsilon);
}

void testSumofTwoExpontentials() {
    std::vector<double> x0 = {0.0, 20.0};
    std::vector<double> x1 = {1.0, 10.0};
    double f0 = two_exponentials(x0);
    double f1 = two_exponentials(x1);

    std::cout << "sumof2exponentials init = " << f0 << std::endl;
    std::cout << "sumof2exponentials min  = " << f0 << std::endl;
    assert(std::abs(f0 - 2.087) < epsilon);
    assert(std::abs(f1) == 0.0);
}

int main() {
    testRosenbrockFunction();
    testWoodsFunction();
    testPowellQuartic();
    testFletcherPowellHelicalValley();
    testFletcherPowellTrigonometric();

    testRosenbrockGradient();
    testWoodsGradient();
    testPowellQuarticGradient();
    testFletcherPowellHelicalValleyDerivative();
    testFletcherPowellTrigonometricDerivative();

    //testThermisterProblem();
    //testSumofTwoExpontentials();

    return 0;
}
