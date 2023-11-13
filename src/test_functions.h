#pragma once
#include "../lib/utility.h"

extern int funev, gradev;

// Rosenbrock and its derivative with regularized variations
double rosenbrock(const std::vector<double>& x);
std::vector<double> rosenbrock_gradient(const std::vector<double>& x);
double rosenbrock_regularized(const std::vector<double>& x);
std::vector<double> rosenbrock_gradient_regularized(const std::vector<double>& x);

// Woods and its derivative
double woods(std::vector<double>& x);
std::vector<double> woods_derivative(const std::vector<double>& x);

// Powell Quartic with Singular Hessian
double powell_quartic(std::vector<double>& x);
std::vector<double> powell_quartic_derivative(const std::vector<double>& x);

// Fletcher and Powell 3 Variable Helical Valley
double helical_valley(std::vector<double>& x);
std::vector<double> helical_valley_derivative(std::vector<double>& x);

// Fletcher Powell Trigonometric Function
double fletcher_powell_trig(std::vector<double>& x0);
std::vector<double> fletcher_powell_trig_derivative(const std::vector<double>& x);

double thermister(std::vector<double>& x);
double two_exponentials(std::vector<double>& x);
double chemical_equilibrium(std::vector<double>& x);
double heat_conduction(std::vector<double>& x);

double randomValue(double lower, double upper);
double rastrigin(std::vector<double>& x);
double ackley(std::vector<double>& x);
double eggholder(std::vector<double>& x);
double goldstein_price(std::vector<double>& x);
