#include "test_functions.h"

int funev = 0, gradev = 0;

// Rosenbrock function
double rosenbrock(const std::vector<double>& x) {
  funev++;
  double sum = 0.0;
  for (size_t i = 0; i < x.size() - 1; ++i) {
      sum += 100 * std::pow(x[i + 1] - std::pow(x[i], 2), 2) + std::pow(1 - x[i], 2);
  }
  return sum;
}

// Gradient of Rosenbrock function
std::vector<double> rosenbrock_gradient(const std::vector<double>& x) {
  gradev++;
  size_t n = x.size();
  std::vector<double> grad(n, 0.0);

  grad[0] = -2*(1-x[0]) - 400 * x[0] * (x[1] - std::pow(x[0], 2));
  for (size_t i = 1; i < n - 1; ++i) {
      grad[i] = -2 * (1 - x[0]) + 200 * (x[i] - std::pow(x[i - 1], 2)) - 400 * x[i] * (x[i + 1] - std::pow(x[i], 2));
  }
  grad[n - 1] = 200 * (x[n - 1] - std::pow(x[n - 2], 2));

  return grad;
}

// Rosenbrock function with L2 Regularization technique to prevent exploding gradient
double rosenbrock_regularized(const std::vector<double>& x) {
  funev++;
  double lambda = 15;
  double sum = 0.0;
  double l2_term = 0.0;
  for (size_t i = 0; i < x.size() - 1; ++i) {
    sum += 100 * std::pow(x[i + 1] - std::pow(x[i], 2), 2) + std::pow(1 - x[i], 2);
    l2_term += std::pow(x[i], 2);
  }
  l2_term += std::pow(x.back(), 2);  // Don't forget the last term
  return sum + lambda * l2_term;
}

std::vector<double> rosenbrock_gradient_regularized(const std::vector<double>& x) {
  gradev++;
  double lambda = 15;
  size_t n = x.size();
  std::vector<double> grad(n, 0.0);

  grad[0] = -400 * x[0] * (x[1] - std::pow(x[0], 2)) - 2 * (1 - x[0]) + 2 * lambda * x[0];
  for (size_t i = 1; i < n - 1; ++i) {
      grad[i] = 200 * (x[i] - std::pow(x[i - 1], 2)) - 400 * x[i] * (x[i + 1] - std::pow(x[i], 2)) - 2 * (1 - x[i]) + 2 * lambda * x[i];
  }
  grad[n - 1] = 200 * (x[n - 1] - std::pow(x[n - 2], 2)) + 2 * lambda * x[n - 1];

  return grad;
}

// Woods Function
double woods(const std::vector<double>& x) {
  funev++;
  double x1 = x[0];
  double x2 = x[1];
  double x3 = x[2];
  double x4 = x[3];
  
  double x1_2 = std::pow(x1, 2);
  double x3_2 = std::pow(x3, 2);

  double term1 = 100*std::pow(x2 - x1_2,2);
  double term2 = std::pow((1 - x1),2) + std::pow((x3 - 1),2);
  double term3 = 90*std::pow((x4 - x3_2),2);
  double term4 = 10.1*(std::pow((x2 - 1),2) + std::pow((x4 - 1),2));
  double term5 = 19.8*(x2 - 1)*(x4 - 1);

  return term1 + term2 + term3 + term4 + term5;
}

// Derivative of Woods function
std::vector<double> woods_derivative(const std::vector<double>& x) {
  gradev++;
  double x1 = x[0];
  double x2 = x[1];
  double x3 = x[2];
  double x4 = x[3];

  double d1 = -400 * x1 * (x2 - std::pow(x1, 2)) - 2 * (1 - x1);
  double d2 = 200 * (x2 - std::pow(x1, 2)) + 20.2 * (x2 - 1) + 19.8 * (x4 - 1);
  double d3 = -360 * x3 * (x4 - std::pow(x3, 2)) + 2 * (x3 - 1);
  double d4 = 180 * (x4 - std::pow(x3, 2)) + 20.2 * (x4 - 1) + 19.8 * (x2 - 1);

  return {d1, d2, d3, d4};
}

// Powell's Quartic Function
double powell_quartic(const std::vector<double>& x) {
  funev++;
  double x1 = x[0];
  double x2 = x[1];
  double x3 = x[2];
  double x4 = x[3];

  double term1 = std::pow(x1 + 10*x2, 2);
  double term2 = 5 * std::pow(x3 - x4, 2);
  double term3 = std::pow(x2 - 2*x3, 4);
  double term4 = 10 * std::pow(x1 - x4, 4);

  return term1 + term2 + term3 + term4;
}

// Powell's Quartic Function
std::vector<double> powell_quartic_derivative(const std::vector<double>& x) {
  gradev++;
  double x1 = x[0];
  double x2 = x[1];
  double x3 = x[2];
  double x4 = x[3];

  // Partial derivatives with respect to x1, x2, x4, and x4
  double term1 = 2 * (x1 + 10 * x2) + 40 * std::pow((x1 - x4),3);
  double term2 = 20 * (x1 + 10 * x2) + 4 * std::pow((x2 - 2 * x3),3);
  double term3 = 10 * (x3 - x4) - 8 * std::pow((x2 - 2 * x3),3);
  double term4 = -10 * (x3 - x4) + 40 * std::pow(x1 - x4,3);

  return {term1, term2, term3, term4};
}


// Fletcher and Powell 3 Variable Helical Valley
double helical_valley(const std::vector<double>& x) {
  funev++;
  double x1 = x[0];
  double x2 = x[1];
  double x3 = x[2];

  const double pi = M_PI;
  double theta;

  // Calculate Theta based on the sign of x1
  /*if (x1 > 0) { // x1 positive
    theta = (1/(2*pi)) * atan2(x2, x1); 
    //std::cout << "x1 > 0, theta = " << theta << std::endl;
} else { // x1 negative
    theta = (1/(2*pi)) + (1/(2*pi)) * atan2(x2, x1);
    //std::cout << "x1 < 0, theta = " << theta << std::endl;
  }*/

  if (x1 < 0) {// negative
    theta = (1/2 * pi) * atan2(x2, x1) + 0.5;
    //std::cout << "x1 > 0, theta = " << theta << std::endl;
  } else { //positive
    theta = (1/2 * pi) * atan2(x2, x1); 
    //std::cout << "x1 < 0, theta = " << theta << std::endl;
  }

  double term1 = 100 * pow(x3 - 10*theta, 2);
  double term2 = pow(sqrt(x1*x1 + x2*x2) - 1, 2);

  return term1 + term2 + x3*x3;
}

// Fletcher and Powell 3 Variable Helical Valley
std::vector<double> helical_valley_derivative(const std::vector<double>& x) {
  gradev++;
  double x1 = x[0];
  double x2 = x[1];
  double x3 = x[2];

  double norm = std::sqrt(x1*x1 + x2*x2);

  // Check to avoid division by zero
  if (norm == 0) {
    return {0.0,0.0,0.0};
  }

  const double pi = M_PI;
  double theta;

  if (x1 < 0) {// negative
    theta = (1/2 * pi) * atan2(x2, x1) + 0.5;
    //std::cout << "x1 > 0, theta = " << theta << std::endl;
  } else { // nonnegative
    theta = (1/2 * pi) * atan2(x2, x1); 
    //std::cout << "x1 < 0, theta = " << theta << std::endl;
  }

  // Partial derivative with respect to x1, x2, and x3
  double term1 =  2 * x1 * (norm - 1) / norm;
  double term2 = 2 * x2 * (norm - 1) / norm;
  double term3 = 200 * (x3 - 10*theta) + 2 * x3;
  return {term1, term2, term3};
}


// Fletcher - Powell Trigonometric function
double fletcher_powell_trig(const std::vector<double>& x0){
  funev++;
  int n = 5 + rand() % 71; // Random n between 5 and 75

  // Initialize x, a, b with random values
  std::vector<double> x(n);
  std::vector<std::vector<double>> a(n, std::vector<double>(n));
  std::vector<std::vector<double>> b(n, std::vector<double>(n));

  for(int i = 0; i < n; i++) {
      // Random value between -pi and pi
      x[i] = -M_PI + (2 * M_PI * (rand() / (double)RAND_MAX));
      for(int j = 0; j < n; j++) {
          // Random values between -100 and 100
          a[i][j] = -100.0 + (200.0 * (rand() / (double)RAND_MAX));
          b[i][j] = -100.0 + (200.0 * (rand() / (double)RAND_MAX));
      }
  }

  n = x.size();
  double sum = 0.0;
  for(int i = 0; i < n; i++) {
      double e_i = 0.0;
      double inner_sum = 0.0;
      for(int j = 0; j < n; j++) {
          double value = a[i][j] * sin(x[j]) + b[i][j] * cos(x[j]);
          e_i += value;
          inner_sum += value;
      }
      sum += (e_i - inner_sum) * (e_i - inner_sum);
  }
  return sum;
}

// Fletcher-Powell Trigonometric function derivative
std::vector<double> fletcher_powell_trig_derivative(const std::vector<double>& x) {
  gradev++;
  int n = x.size();

  // Initialize a, b with random values
  std::vector<std::vector<double>> a(n, std::vector<double>(n));
  std::vector<std::vector<double>> b(n, std::vector<double>(n));

  for(int i = 0; i < n; i++) {
      for(int j = 0; j < n; j++) {
          a[i][j] = -100.0 + (200.0 * (rand() / (double)RAND_MAX));
          b[i][j] = -100.0 + (200.0 * (rand() / (double)RAND_MAX));
      }
  }

  // Gradient vector
  std::vector<double> grad(n, 0.0);

  for(int k = 0; k < n; k++) {
      double gradient = 0.0;
      for(int i = 0; i < n; i++) {
          double e_i = 0.0;
          for(int j = 0; j < n; j++) {
              e_i += a[i][j] * sin(x[j]) + b[i][j] * cos(x[j]);
          }

          double derivative_e_i = a[i][k] * cos(x[k]) - b[i][k] * sin(x[k]);
          gradient += 2 * (e_i - e_i) * derivative_e_i;
      }
      grad[k] = gradient;
  }
  return grad;
}

double randomValue(double lower, double upper) {
    return lower + (upper - lower) * (rand() / (double)RAND_MAX);
}

// Function to calculate the predicted Y_i
double model(double x1, double x2, double x3, double Ti) {
    return x1 * std::exp(x2 / (Ti + x3));
}


double thermister(const std::vector<double>& x) {
  funev++;
  const int n = 16;
  
  // Initialize y_hat, T with random values ???
  std::vector<double> y_hat(n);
  std::vector<double> T(n);
  for(int i = 0; i < n; i++) {
      y_hat[i] = randomValue(0.0, 1.0);
      T[i] = randomValue(0.0, 1.0);
  }

  double x1 = x[0];
  double x2 = x[1];
  double x3 = x[2];
  
  double sum = 0.0;

  for(int i = 0; i < n; i++) {
      double y_i = x1 * exp(x2 / (T[i] + x3));
      sum += (y_i - y_hat[i]) * (y_i - y_hat[i]);
  }

  return sum;
}

// Sum of Two Exponentials
double two_exponentials(const std::vector<double>& x) {
  funev++;
  double x1 = x[0]; 
  double x2 = x[1];

  double sum = 0.0;

  for(int i=1; i<=10; i++) {
    double ti = 0.1 * i;
    
    double term1 = std::exp(-x1*ti);
    double term2 = std::exp(-x2*ti);
    double term3 = std::exp(ti) - std::exp(-10*ti);
    double subtotal = std::pow((term1 - term2) - term3, 2);
    sum += subtotal;
  }
  return sum;
}

// Chemical Equilibrium Problem 
double chemical_equilibrium(const std::vector<double>& x) {
  funev++;
  double x1 = x[0];
  double x2 = x[1];
  double x3 = x[2];

  double term1 = pow((1 - x1 - x2)*(1 - x3 - x1) - 4*pow(x1,2)/549000, 2); 
  double term2 = pow((1 - x1 - x2)*(1 - x2 - x3) - 4*pow(x2,2)/362, 2);
  double term3 = pow((1 - x2 - x3)*(1 - x3 - x1) - 4*pow(x3,2)/3.28, 2);

  return term1 + term2 + term3;
}

// Heat Conduction Problem 
double heat_conduction(const std::vector<double>& x ) {
  funev++;
  double x1 = x[0];
  double x2 = x[1]; 
  double x3 = x[2];
  double x4 = x[3];

  double term1 = pow(2*(x2 + x3 - 4*x1) + 20 - 1.5*x1 + pow(x1,2)/20, 2);
  double term2 = pow(2*(x1 - 3*x3 + x4) + 20 - 1.5*x3 + pow(x3,2)/20, 2);
  double term3 = pow(2*(x2 - x4) + 20 - 1.5*x2 + pow(x2,2)/20, 2);
  double term4 = pow(2*(x1 + x3 - 2*x4) + 20 - 1.5*x4 + pow(x4,2)/20, 2);

  return term1 + term2 + term3 + term4;

}

// simulated annealing
// 

double rastrigin(const std::vector<double>& x) {
  funev++;
  double sum = 0;

  for (int i = 0; i < x.size(); i++) {
    sum += x[i] * x[i] - 10 * cos(2 * M_PI * x[i]); 
  }
  return 10 * x.size() + sum;
}

double ackley(const std::vector<double>& x) {
  funev++;
    double sum1 = 0;
    double sum2 = 0;
    for (int i = 0; i < x.size(); i++) {
        sum1 += x[i] * x[i];
        sum2 += cos(2 * M_PI * x[i]);
    }
    return -20 * exp(-0.2 * sqrt(sum1 / x.size())) - exp(sum2 / x.size()) + 20 + M_E; 
}

double eggholder(const std::vector<double>& x) {
  funev++;
  return -x[1] * sin(sqrt(abs(x[0] + x[1] + 47))) 
         - x[0] * sin(sqrt(abs(x[0] - (x[1] + 47))));
}

double goldstein_price(const std::vector<double>& x) {
  funev++;
    double x1 = x[0];
    double x2 = x[1];

    double term1 = (1 + pow(x1 + x2 + 1, 2) * (19 - 14*x1 + 3*pow(x1, 2) - 14*x2 + 6*x1*x2 + 3*pow(x2, 2)));
    double term2 = (30 + pow(2*x1 - 3*x2, 2) * (18 - 32*x1 + 12*pow(x1, 2) + 48*x2 - 36*x1*x2 + 27*pow(x2, 2)));

    return term1 * term2;
}
