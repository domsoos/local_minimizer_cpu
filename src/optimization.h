#pragma once
#include <deque>
#include "utility.h"
#include "test_functions.h"
#include "genetic.h"


// Line Search methods
double simple_backtracking(std::function<double(std::vector<double> &)> func, std::vector<double> x, std::vector<double> p, double alpha, double tau);
double line_search_simple(std::function<double(std::vector<double> &)> func, std::vector<double> x, std::vector<double> p, double alpha, double tau);
double line_search(std::function<double(std::vector<double> &)> func, std::vector<double> x, std::vector<double> p, double alpha, double c, double tau);
double cubicInterpolationLineSearch(std::function<double(std::vector<double> &)> func, std::vector<double> x, std::vector<double> s, double f0);
double mnLineSearch(std::function<double(std::vector<double> &)> func, std::vector<double> x, std::vector<double> p, double gdel);
double quadratic_line_search(std::function<double(std::vector<double> &)> func, std::vector<double> x, std::vector<double> p, double alpha, double c, double tau, int maxiter);
double safe_divide(double numerator, double denominator, double default_value);
double armijoCurvature(std::function<double(std::vector<double> &)> func, std::vector<double> x, std::vector<double> p, std::vector<double> grad, double alpha, double tau);

double dlib_line_search(std::function<double(std::vector<double> &)> func,double f0,double d0,double rho,double sigma,double min_f,const int max_iter,std::vector<double> &x,std::vector<double> &p, double lambda);
void printVectors(int i, double fnew,std::vector<double> x,std::vector<double> delta_x,std::vector<double> delta_g,double alpha);

double directional_derivative(const std::vector<double> &grad, const std::vector<double> &p);
double backtracking_line_search(std::function<double(std::vector<double>&)> func,const std::vector<double>& x,const std::vector<double>& g,const std::vector<double>& p,double alpha,double beta,double c, double step_size_lower_limit,double default_alpha);
double wolfe_line_search(std::function<double(const std::vector<double>&)> func,std::function<std::vector<double>(const std::vector<double>&)> grad_func,const std::vector<double>& x,const std::vector<double>& gradient,std::vector<double>& direction,double alpha = 1.0,double c1=1e-4,double c2=0.9);
std::vector<double> simulated_annealing(std::function<double(std::vector<double> &)> func, std::vector<double> initial_state, double initial_temp, double alpha, double temp_threshold, int max_iter);

std::pair<std::vector<double>, std::vector<double>> lbfgsb_step(std::function<double(std::vector<double> &)> func,double alpha, std::vector<double> p,std::string algorithm,
    const std::vector<double>& x, const std::vector<double>& g,
    const double lower,const double upper, // bounds.first = lower else upper
    std::deque<std::vector<double>>& s_history, std::deque<std::vector<double>>& y_history, // histories
    std::deque<double>& rho_history, double& gamma_k, int m); // m = history size, gamma_k = scaling factor
std::vector<double> lbfgsb_update(const std::vector<double> &g, std::deque<std::vector<double>> &s_history, 
    std::deque<std::vector<double>> &y_history, std::deque<double> &rho_history,
    double gamma_k, int m);

void bfgs_update(std::vector<std::vector<double>>& H, std::vector<double> delta_x, std::vector<double> delta_g, double delta_dot);
void dfp_update(std::vector<std::vector<double>>& H, std::vector<double> delta_x, std::vector<double> delta_g);

double optimize(std::function<double(std::vector<double> &)> func, std::vector<double> x0, std::string algorithm,const double tol,const int max_iter,const double lower, double upper);
double minimize(std::function<double(std::vector<double> &)> func, std::vector<double> x0, std::string name,const int pop_size,const int max_gens,const int dim, std::string algorithm, double lower, double upper);