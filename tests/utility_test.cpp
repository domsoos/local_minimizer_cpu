#include "utility.h"
#include <cassert>
#include <cmath>

void test_rand_double() {
    double min = -5.0;
    double max = 5.0;
    double val = rand(min, max);
    std::cout << "rand_double() = "<< val << std::endl;
    assert(val >= min && val <= max);
}

void test_rand_int() {
    int min = -5;
    int max = 5;
    int val = rand(min, max);
    std::cout << "rand_int() = " << val << std::endl;
    assert(val >= min && val <= max);
}

void test_measure_memory() {
    long memory = measure_memory();
    std::cout << "measure_memory() = " << memory << std::endl;
    assert(memory > 0); // Assuming the function will always return a positive value
}

void test_scale_vector() {
    std::vector<double> v{1.0, 2.0, 3.0};
    double scalar = 2.0;
    std::vector<double> result = scale_vector(v, scalar);
    assert(result[0] == 2.0 && result[1] == 4.0 && result[2] == 6.0);
}

void test_add_vectors() {
    std::vector<double> v1{1.0, 2.0, 3.0};
    std::vector<double> v2{4.0, 5.0, 6.0};
    std::vector<double> result = add_vectors(v1, v2);
    assert(result[0] == 5.0 && result[1] == 7.0 && result[2] == 9.0);
}

int main() {
    test_rand_double();
    test_rand_int();
    test_measure_memory();
    test_scale_vector();
    test_add_vectors();

    return 0;
}
