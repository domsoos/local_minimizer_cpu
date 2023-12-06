#include "utility.h"
#include "test_functions.h"
#include "optimization.h"

#include <sys/resource.h>
#include <dlib/optimization.h>

#include <fstream>
#include <sstream>

//int funev = 0, gradev = 0;

using FunctionType = std::function<double(const std::vector<double>&)>;
using DerivativeType = std::function<std::vector<double>(const std::vector<double>&)>;

typedef dlib::matrix<double,2,1> column_vector2D;
typedef dlib::matrix<double,3,1> column_vector3D;
typedef dlib::matrix<double,4,1> column_vector4D;
typedef dlib::matrix<double,30,1> column_vector30D;

void printdlibvector(dlib::matrix<double> c) {
    for(int i=0;i<c.size();i++) {
        std::cout << c(i) << " ";
    }
    std::cout << std::endl;
}

double rosenbrock_dlib_wrapper(const column_vector2D& m)
{
  funev++;
  return rosenbrock({m(0), m(1)});
}

double 
woods_dlib_wrapper(const column_vector4D& m)
{
    funev++;
    return woods({m(0),m(1),m(2),m(3)});
}

double
powell_dlib_wrapper(const column_vector4D& m)
{
    funev++;
    return powell_quartic({m(0),m(1),m(2),m(3)});
}

double 
helical_dlib_wrapper(const column_vector3D& m)
{
    funev++;
    return helical_valley({m(0),m(1),m(2)});
}

double 
trig_dlib_wrapper(const column_vector30D& m)
{
    funev++;
    std::vector<double> v(30);
    for(int i=0; i<30;i++) {
        v.push_back(m(i));
    }
    return fletcher_powell_trig(v);
}

struct  ds {
    double min;
    long it;
};

template <typename MatrixType>
std::vector<double> dlib_optimization(MatrixType x, std::function<double(const MatrixType&)> fun_wrapper) {
    std::cout << "fun init = "<< fun_wrapper(x) << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    dlib::find_min_using_approximate_derivatives(dlib::bfgs_search_strategy(),
                                    dlib::objective_delta_stop_strategy(1e-6),
                                    fun_wrapper, x, 0.0);
    //std::cout << "fmin = " << belo[0] << "\nit = " << belo[1] << std::endl;
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    long time = duration.count();
    printdlibvector(x);
    std::cout << "min = " << std::abs(fun_wrapper(x)) << std::endl;
    std::cout << "\ntime = " << time << " ms" << std::endl;
    std::cout << "funev =  " << funev << std::endl;
    return {static_cast<double>(time), fun_wrapper(x)};
}

int main() {
    int dim;
    //long before, after, result;
    double error;
    char going;
    bool done = false;
    funev = 0, gradev = 0;
    double lower, upper;
    double x = distr(eng);
    std::vector<double> trig_vec(30);
    std::vector<double> true_trig(30);
    for(int i=0;i<30;i++){
        trig_vec.push_back(distr(eng));
        true_trig.push_back(0.0);
    }

    std::vector<std::vector<double>> initialPoints = {
        {-1.2,1.0}, // Rosenbrock 
        {-3.0, -1.0, -3.0, -1.0}, // Woods 
        {3.0, -1.0, 0.0, 1.0}, // Powell
        //{-1.0, 0.0, 0.0}, // FP Helical Valley
        {-6.88376800380764878e+05,3.00693715986630297e+05,5.11036303162183205e+05},
        //{0.362191,10.932092, 10.0}, //-3.82465e-07,-3.90997e-07    1.3746e-10
        trig_vec // FP Trig
    };

    std::vector<FunctionType> functions = {
        rosenbrock,
        woods,
        powell_quartic,
        helical_valley,
        fletcher_powell_trig,
    };

    std::vector<DerivativeType> derivatives = {
        rosenbrock_gradient,
        woods_derivative,
        powell_quartic_derivative,
        helical_valley_derivative,
        fletcher_powell_trig_derivative
    };

    std::vector<std::string> const names = {
        "Rosenbrock",
        "Woods",
        "Powell Quartic",
        "Fletcher Powell Helical Valley",
        "Fletcher Powell Trigonometric"
    };

    std::vector<std::vector<double>> const trueMin = {// last in vec is funmin
        {1.0,1.0,0.0}, // Rosen
        {1.0,1.0,1.0,1.0,0.0}, // Woods
        {0.0,0.0,0.0,0.0,0.0}, // Powell
        {1.0,0.0,0.0,0.0}, // Helical
        true_trig, // Trig
    };

    std::string algorithm;
    while(!done) {
        algorithm = getAlgorithm();
        for(int it=0;it<functions.size()-1;it++) {
            funev = 0, gradev = 0;
            if (algorithm  == "lbfgsb"){
                std::cout<<"Enter lower bound: ";
                std::cin>>lower;
                std::cout<<"Enter upper bound: ";
                std::cin>>upper;
            }// end if lbfgsb

            std::vector<double> x0 = initialPoints[it];
            std::function<double(const std::vector<double>&)> fun = functions[it];
            std::function<std::vector<double>(const std::vector<double>&)> der = derivatives[it];

            std::cout << "\n\n\t" << names[it] <<" Function Optimization\n\n\n\t\t~~~Celsius Optimization~~~\n";
            
            std::cout << "x0 = ";
            for(int i=0;i<x0.size();i++) {std::cout << x0[i] << " ";}
            std::cout << "\ninit fun = " << fun(x0) << std::endl;

            auto start = std::chrono::high_resolution_clock::now();
            error = minimize(fun,der,x0,names[it],0,0,dim, algorithm, lower, upper);
            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            long time = duration.count();
            std::cout << "\ntime: " << time << " ms" << std::endl;

            std::cout << "\nGlobal Minimum(";
            for(int i=0;i<trueMin[it].size()-1;i++) {
                std::cout << trueMin[it][i] << ",";
            }
            std::cout << ") = " << trueMin[it][-1] << "\nError = " << error << std::endl;
            std::cout << "Function Evaluations: " << funev << "\tGradient Evaluations: " << gradev << std::endl << std::endl;
        
            funev = 0, gradev = 0;
            std::cout << "\n\t\t~~~DLIB OPTIMIZATION~~~" << std::endl;
            if(names[it] == "Rosenbrock") {
                column_vector2D r0 = {-1.2,1.0};
                std::function<double(const column_vector2D&)> func_wrapper = rosenbrock_dlib_wrapper;
                std::vector<double> belo = dlib_optimization(r0,func_wrapper);
            }
            else if(names[it] ==  "Woods") {
                column_vector4D w0 = {-3.0, -1.0, -3.0, -1.0};
                std::function<double(const column_vector4D&)> func_wrapper = woods_dlib_wrapper;
                std::vector<double> belo = dlib_optimization(w0,func_wrapper);
            }
            else if(names[it] == "Powell Quartic") {
                column_vector4D p0 = {3.0, -1.0, 0.0, 1.0};
                std::function<double(const column_vector4D&)> func_wrapper = powell_dlib_wrapper;
                std::vector<double> belo = dlib_optimization(p0,func_wrapper);
            }
            else if(names[it] == "Fletcher Powell Helical Valley") {
                column_vector3D h0 =  {-6.88376800380764878e+05,3.00693715986630297e+05,5.11036303162183205e+05};
                std::function<double(const column_vector3D&)> func_wrapper = helical_dlib_wrapper;
                std::vector<double> belo = dlib_optimization(h0,func_wrapper);         
            }
            else if(names[it] == "Fletcher Powell Trigonometric") {
                column_vector30D t0 = {distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng),distr(eng)};
                std::function<double(const column_vector30D&)> func_wrapper = trig_dlib_wrapper;
                std::vector<double> belo = dlib_optimization(t0,func_wrapper);  
            }
            // Save belo to the file belo[0] = time, belo[1] =
        }
        done = true;
        /*
        std::cout << "\n\nKeep going? ";
        std::cin >> going;
        if(std::tolower(going) != 'y'){
            std::cout << "Goodbye!"<<std::endl;
            done = true;
        }
        */
    }//end while
    return 0;
}