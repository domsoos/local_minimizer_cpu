#pragma once
#include "../lib/utility.h"

extern int funev, gradev;

// Rosenbrock and its derivative with regularized variations
double rosenbrock(const std::vector<double>& x);
std::vector<double> rosenbrock_gradient(const std::vector<double>& x);
double rosenbrock_regularized(const std::vector<double>& x);
std::vector<double> rosenbrock_gradient_regularized(const std::vector<double>& x);

// Woods and its derivative
double woods(const std::vector<double>& x);
std::vector<double> woods_derivative(const std::vector<double>& x);

// Powell Quartic with Singular Hessian
double powell_quartic(const std::vector<double>& x);
std::vector<double> powell_quartic_derivative(const std::vector<double>& x);

// Fletcher and Powell 3 Variable Helical Valley
double helical_valley(const std::vector<double>& x);
std::vector<double> helical_valley_derivative(const std::vector<double>& x);

// Fletcher Powell Trigonometric Function
double fletcher_powell_trig(const std::vector<double>& x0);
std::vector<double> fletcher_powell_trig_derivative(const std::vector<double>& x);

double thermister(const std::vector<double>& x);
double two_exponentials(const std::vector<double>& x);
double chemical_equilibrium(const std::vector<double>& x);
double heat_conduction(const std::vector<double>& x);

double randomValue(double lower, double upper);
double rastrigin(const std::vector<double>& x);
double ackley(const std::vector<double>& x);
double eggholder(const std::vector<double>& x);
double goldstein_price(const std::vector<double>& x);
