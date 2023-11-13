#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <random>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cmath>
#include <limits>

double rand(double min, double max);
int rand(int min, int max);
long measure_memory();

extern double global_min;
extern std::vector<double> best_params;

extern std::mt19937 eng;
extern std::uniform_real_distribution<> distr;

void print_vector(std::vector<double> v);
void printVectors(int i, double fnew,std::vector<double> x,std::vector<double> delta_x,std::vector<double> delta_g,double alpha);

bool isValidDouble(double val);
std::vector<double> scale_vector(const std::vector<double> &v1, double scalar);
std::vector<double> add_vectors(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> subtract_vectors(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> matvec_product(const std::vector<std::vector<double>>& m, const std::vector<double>& v);
double dot_product(const std::vector<double>& v1, const std::vector<double>& v2);
double norm(const std::vector<double>& v);
std::vector<std::vector<double>> matmul(const std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2);
std::vector<std::vector<double>> outer_product(const std::vector<double>& v1, const std::vector<double>& v2);
std::vector<double> gradient(std::function<double(std::vector<double> &)> func, std::vector<double> x, double h);

std::vector<double> rosenbrock_multi_gradient(const std::vector<double>& x);