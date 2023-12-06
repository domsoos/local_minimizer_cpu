#include "utility.h"

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution(-5.0, 5.0);

std::random_device rd;
std::mt19937 eng(rd());
std::uniform_real_distribution<> distr(-M_PI, M_PI);


double rand(double min, double max) {
  return distribution(generator); 
}

int rand(int min, int max) {
  std::uniform_int_distribution<int> dist(min, max);
  return dist(generator);
}

long measure_memory() {
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    return r_usage.ru_maxrss;
}

double global_min = std::numeric_limits<double>::max();
std::vector<double> best_params;


void print_vector(std::vector<double> v) {
    for(int i=0;i<v.size();i++) {
        std::cout << "x["<<i<<"]: " << v[i] << std::endl;
    }
}
void printVectors(int i, double fnew,std::vector<double> x,std::vector<double> delta_x,std::vector<double> delta_g,double alpha){
    std::cout<< "it #" << i << std::endl << "fun = " << fnew<<std::endl;
    print_vector(x);
    std::cout << "delta_x: ";
    print_vector(delta_x);
    std::cout << "\ndelta_g: ";
    print_vector(delta_g);
    std::cout << "\n" << std::endl;
    std::cout << "current alpha =  " << alpha <<std::endl;
}

std::string getAlgorithm() {
    char answer, lbfgs, box;
    std::cout << "\nChoose the algorithm"<<std::endl;;
    std::cout << "Use BFGS? ";
    std::cin >> answer;
    if (std::tolower(answer) == 'y') {
        std::cout << "Limited memory variant - L-BFGS? ";
        std::cin >> lbfgs;
        if (std::tolower(lbfgs) == 'y') {
            std::cout << "bounded? ";
            std::cin >> box;
            if(std::tolower(box) == 'y') {
                return "lbfgsb";
            } else {
                return "lbfgs";
            }
        } else {
            return "bfgs";
        }
    }
    else { return "dfp";}
}

bool isValidDouble(double val) {
    if (std::isnan(val) || std::isinf(val)) {std::cout << "invalid value: "<<val<<std::endl;return false; }
    else {return true;}
} 

/*** Utility Functions for Matrix-Vector-Scalar operations  ***/
// Scale a vector
std::vector<double> scale_vector(const std::vector<double> &v1, double scalar) {
    std::vector<double> result(v1.size());
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] * scalar;
    }//end for
    return result;
}
// Vector addition
std::vector<double> add_vectors(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vector sizes do not match for addition");
    }// end if

    std::vector<double> result(v1.size(), 0.0);
    for (size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] + v2[i];
    }// end for
    return result;
}// end add_vectors 

// Vector subtraction
std::vector<double> subtract_vectors(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size()) {
        throw std::invalid_argument("Vector sizes do not match for subtraction");
    }// end if

    std::vector<double> result(v1.size());
    for (std::size_t i = 0; i < v1.size(); ++i) {
        result[i] = v1[i] - v2[i];
    }// end for

    return result;
}// end subtract_vectors

// Matrix-vector multiplication
std::vector<double> matvec_product(const std::vector<std::vector<double>>& m, const std::vector<double>& v) {
    int rows = m.size();
    int cols = m[0].size();

    if (cols != v.size()) {
        throw std::runtime_error("Invalid dimensions for matrix-vector multiplication");
    }// end if

    std::vector<double> result(rows, 0.0);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i] += m[i][j] * v[j];
        }// end inner for
    }// end outer for
    return result;
} // end matvec_product

// Dot product of two vectors
double dot_product(const std::vector<double>& v1, const std::vector<double>& v2) {
    double result = 0.0;
    for (int i = 0; i < v1.size(); i++) {
        result += v1[i] * v2[i];
    }// end for
    return result;
}// end dot_product

// Norm of a vector
double norm(const std::vector<double>& v) {
    return std::sqrt(dot_product(v, v));
}// end norm

// Matrix multiplication
std::vector<std::vector<double>> matmul(const std::vector<std::vector<double>>& m1, const std::vector<std::vector<double>>& m2) {
    int rows1 = m1.size();
    int cols1 = m1[0].size();
    int rows2 = m2.size();
    int cols2 = m2[0].size();

    if (cols1 != rows2) {
        throw std::runtime_error("Invalid dimensions for matrix multiplication");
    }// end if

    std::vector<std::vector<double>> result(rows1, std::vector<double>(cols2, 0));
    // 6 x 6 matrix multiplication
    for (int i = 0; i < rows1; i++) { // for each row in the matrix
        for (int j = 0; j < cols2; j++) { // for each column in the 2nd matrix
            // For each row in the first matrix and column in the second matrix, calculate the dot product of the corresponding row and column
            // The dot product is calculated by multiplying the corresponding elements together and then summing the result
            for (int k = 0; k < cols1; k++) {
                result[i][j] += m1[i][k] * m2[k][j];
            }// end innnner
        }// end inner
    }// end outer
    return result;
}// end matmul

/* Outer product of two vectors
 * The outer product of two vectors is a matrix where each element is the product of an element from the first vector and an element from the second vector.
 */
std::vector<std::vector<double>> outer_product(const std::vector<double>& v1, const std::vector<double>& v2) {
    int size1 = v1.size();
    int size2 = v2.size();

    std::vector<std::vector<double>> result(size1, std::vector<double>(size2, 0));
    for (int i = 0; i < size1; i++) {
        for (int j = 0; j < size2; j++) {
            result[i][j] = v1[i] * v2[j];
        }// end inner
    }// end outer
    return result;
}// end outer_product

// Finite difference gradient calculation
std::vector<double> gradient(std::function<double(std::vector<double> &)> func, std::vector<double> x, double h) {
    std::vector<double> grad(x.size(), 0);
    for (int i = 0; i < x.size(); i++) {
        double temp = x[i];
        x[i] = temp + h;
        double fp = func(x);
        x[i] = temp - h;
        double fm = func(x);
        grad[i] = (fp - fm) / (2 * h);
        x[i] = temp;
    }// end for
    return grad;
}// end gradient

std::vector<double> rosenbrock_multi_gradient(const std::vector<double>& x) {
    int n = x.size();
    std::vector<double> grad(n, 0);

    // Special case for the first element
    grad[0] = -400 * x[0] * (x[1] - x[0] * x[0]) + 2 * (1 - x[0]);

    // Special case for the last element
    grad[n - 1] = 200 * (x[n - 1] - x[n - 2] * x[n - 2]);

    // For the elements in between
    for (int i = 1; i < n - 1; ++i) {
        grad[i] = -400 * x[i] * (x[i + 1] - x[i] * x[i]) + 200 * (x[i] - x[i - 1] * x[i - 1]) + 2 * (1 - x[i]);
    }

    return grad;
}