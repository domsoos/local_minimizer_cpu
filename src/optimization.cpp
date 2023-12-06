#include "optimization.h"

/* Line Minimization from the original paper 
 * that use cubic inteprolation
 */
double cubicInterpolationLineSearch(std::function<double(std::vector<double> &)> func, std::vector<double> x, std::vector<double> s, double f0) {
    double f = func(x);
    std::vector<double> g = gradient(func, x, 1e-7);
    double q = dot_product(g, s);

    // calculate lambda
    double lambda = std::min(1.0, std::abs(-2 * (f - f0) / q));
    std::vector<double> x_prime = x;
    for (auto& val : x_prime) {
        val += lambda;
    }// end for

    double f_prime = func(x_prime);
    std::vector<double> g_prime = gradient(func, x_prime, 1e-7);
    double q_prime = dot_product(g_prime, s);

    double z = 3 * (f - f_prime) / lambda + q + q_prime;
    double a =  z * z - q * q_prime;
    double w;
    if (a > 0) {
        w = std::sqrt(z * z - q * q_prime);    
    } else {
        std::cout << "can't take square root"<< std::endl;
        w = 0;
    }
    //double w = std::sqrt(z * z - q * q_prime);
    //std::cout <<"\n\nlambda = " << lambda << std::endl;
    //std::cout << "q = " << q << "\nq'= "<<q_prime << "\nz * z - q * q_prime = " << a;
    if (z * z < q * q_prime){// || q > 0 || q_prime < 0) {
        throw std::invalid_argument("Conditions not met for line minimization.");
    }// end if

    double alpha = (z + w - q) / (q_prime + 2 * w - q);
    //std::cout << "\nalpha = " << alpha;
    return alpha;
}// end cubicInterpolationLineSearch

double simple_backtracking(std::function<double(std::vector<double> &)> func, std::vector<double> x, std::vector<double> p, double alpha, double tau) {
    double fx = func(x);
    std::vector<double> x_alpha_p = x;  // new proposed point in the parameter space

    // Adjust x_alpha_p along the search direction p by the initial step size of alpha
    for (int i = 0; i < x.size(); i++) {
        x_alpha_p[i] += alpha * p[i];
    }// end for

    // If the function value at the new proposed point is not less than at the current point,
    //    alpha is reduced by a factor of tau to backtrack along the search direction p
    while (func(x_alpha_p) >= fx) {
        alpha *= tau;

        // Adjust x_alpha_p to the new value of alpha
        for (int i = 0; i < x.size(); i++) {
            x_alpha_p[i] = x[i] + alpha * p[i];
        }// end for
    }// end while
    return alpha;
}

double quadratic_line_search(std::function<double(std::vector<double> &)> func,
                             std::vector<double> x,
                             std::vector<double> p,
                             double alpha = 1.0,
                             double c = 0.5,
                             double tau = 0.5,
                             int maxiter = 12) {
    // Constants for the search
    const double SLAMBG = 5.;
    //const double ALPHA = 2.;
    const double TOLER = 0.05;

    int niter = 0;
    double slam = 1.;
    double toler8 = TOLER;
    double slamax = SLAMBG;
    double fvmin = func(x);
    double xvmin = 0.;

    while (niter < maxiter) {
        std::vector<double> new_point = x;
        for (size_t i = 0; i < x.size(); i++) {
            new_point[i] += slam * p[i];
        }
        //std::vector<double> g = gradient(func, x, 1e-7);
        std::vector<double> g = rosenbrock_multi_gradient(x);

        double f_new = func(new_point);

        if (f_new < fvmin) {
            fvmin = f_new;
            xvmin = slam;
        }

        // Check the Armijo condition
        if (f_new > fvmin + c * slam * dot_product(g, p)) {
            slam *= tau;
        } else {
            // If it's not met, use the quadratic model to update slam
            double a = (f_new - fvmin) / (slam * slam);
            slam = -dot_product(g, p) / (2 * a);
        }

        // Enforce step size bounds
        if (slam < toler8) slam = toler8;
        if (slam > slamax) slam = slamax;

        niter++;
    }

    return xvmin;
}

double safe_divide(double numerator, double denominator, double default_value = std::numeric_limits<double>::max()) {
    if (std::abs(denominator) < 1e-10) {
        return default_value;
    }
    return numerator / denominator;
}


/* Armijo Backtracking Line Search to find an appropriate step size alpha 
 * along the direction p at point x for the function func.
 */
double armijoCurvature(std::function<double(std::vector<double> &)> func, std::vector<double> x, std::vector<double> p, std::vector<double> grad, double alpha, double tau) {
    double c1 = 1e-3;
    double c2 = 0.1;
    double fx = func(x);
    std::vector<double> x_alpha_p = x;  // new proposed point in the parameter space

    // Adjust x_alpha_p along the search direction p by the initial step size of alpha
    for (int i = 0; i < x.size(); i++) {
        x_alpha_p[i] += alpha * p[i];
    }// end for

    // While the Armijo condition or the Curvature condition is not met,
    //    alpha is reduced by a factor of tau to backtrack along the search direction p
    while (func(x_alpha_p) > fx + c1 * alpha * dot_product(grad, p) || dot_product(grad, p) < c2 * dot_product(grad, p)) {
        alpha *= tau;

        // Adjust x_alpha_p to the new value of alpha
        for (int i = 0; i < x.size(); i++) {
            x_alpha_p[i] = x[i] + alpha * p[i];
        }// end for
    }// end while
    return alpha;
}// end armijoCurvature

//alpha = armijoCurvature(func, x, p, g, alpha, 0.2);
//alpha = wolfe_line_search(func, rosenbrock_gradient_regularized, x, g, p);//,1.0, 1e-4,0.9);
// Wolf Condition involves both the armijo condition and the curvature condition
double wolfe_line_search(std::function<double(std::vector<double> &)> func,std::function<std::vector<double>(std::vector<double>&)> grad_func,const std::vector<double>& x,const std::vector<double>& gradient,std::vector<double>& direction,double alpha = 1.0,double c1 = 1e-4,double c2 = 0.9) {
    int maxiter = 100;
    int i = 0;
    std::vector<double> x0 = x;
    double f_x = func(x0);
    double grad_dot_dir = dot_product(gradient, direction);

    while (true) {
        i++;
        if(i >= maxiter) {break;}
        std::vector<double> x_alpha = add_vectors(x,scale_vector(direction, alpha));
        double f_x_alpha = func(x_alpha);
        if (f_x_alpha > f_x + c1 * alpha * grad_dot_dir) {
            // Armijo condition not satisfied
            alpha *= 0.5;
            //std::cout << "still in armijo";
            continue;
        }

        std::vector<double> grad_x_alpha = grad_func(x_alpha);
        double wolfe_condition = dot_product(grad_x_alpha, direction);
        if (wolfe_condition < c2 * grad_dot_dir) {
            // Curvature condition not satisfied
            //std::cout <<"made it to curvature";
            alpha *= 2.0;
            continue;
        }
        break;
    }
    return alpha;
}

double directional_derivative(const std::vector<double> &grad, const std::vector<double> &p) {
    double d = 0.0;
    for (size_t i = 0; i < grad.size(); ++i) {
        d += grad[i] * p[i];
    }
    return d;
}

template <typename T>
T put_in_range(const T& a, const T& b, const T& val) {
    if (a < b) {
        if (val < a) return a;
        else if (val > b) return b;
    } else {
        if (val < b) return b;
        else if (val > a) return a;
    }
    return val;
}

inline double put_in_range(const double& a, const double& b, const double& val)
{
    return put_in_range<double>(a, b, val);
}

inline double poly_min_extrap(double f0, double d0, double f1, double d1, double limit = 1) {
    const double n = 3 * (f1 - f0) - 2 * d0 - d1;
    const double e = d0 + d1 - 2 * (f1 - f0);

    double temp = std::max(n * n - 3 * e * d0, 0.0);

    if (temp < 0) return 0.5;

    temp = std::sqrt(temp);

    if (std::abs(e) <= std::numeric_limits<double>::epsilon()) return 0.5;

    double x1 = (temp - n) / (3 * e);
    double x2 = -(temp + n) / (3 * e);

    double y1 = f0 + d0 * x1 + n * x1 * x1 + e * x1 * x1 * x1;
    double y2 = f0 + d0 * x2 + n * x2 * x2 + e * x2 * x2 * x2;

    double x = (y1 < y2) ? x1 : x2;

    return put_in_range(0.0, limit, x);
}

inline double poly_min_extrap(double f0, double d0, double f1) {
    const double temp = 2 * (f1 - f0 - d0);
    if (std::abs(temp) <= d0 * std::numeric_limits<double>::epsilon()) return 0.5;

    const double alpha = -d0 / temp;
    return put_in_range(0.0, 1.0, alpha);
}

inline double poly_min_extrap(double f0, double d0, double x1, double f_x1, double x2, double f_x2) {
    const double aa2 = x2 * x2;
    const double aa1 = x1 * x1;
    double temp = aa2 * aa1 * (x1 - x2);

    if (temp == 0 || std::fpclassify(temp) == FP_SUBNORMAL) return x1 / 2.0;

    double m11 = aa2, m12 = -aa1;
    double m21 = -aa2 * x2, m22 = aa1 * x1;
    double v1 = f_x1 - f0 - d0 * x1;
    double v2 = f_x2 - f0 - d0 * x2;

    double a = (m11 * v1 + m12 * v2) / temp;
    double b = (m21 * v1 + m22 * v2) / temp;

    temp = b * b - 3 * a * d0;

    if (temp < 0 || a == 0) return (f0 < f_x2) ? 0 : x2;

    return put_in_range(0.0, x2, (-b + std::sqrt(temp)) / (3 * a));
}


double dlib_line_search(std::function<double(std::vector<double> &)> f, std::function<std::vector<double>(const std::vector<double>&)> der, double f0,double d0,double rho,double sigma,double min_f,const int max_iter,std::vector<double> &x,std::vector<double> &p, double lambda)
    {
    //std::function<std::vector<double>(const std::vector<double>&)> der,
    if (!(0 < rho && rho < sigma && sigma < 1 && max_iter > 0)) {
        std::cerr << "Invalid arguments provided to line_search" << std::endl;
        throw std::invalid_argument("Invalid arguments provided to line_search");
        return -1;
    }

    // The bracketing phase of this function is implemented according to block 2.6.2 from
    // the book Practical Methods of Optimization by R. Fletcher. The sectioning 
    // phase is an implementation of 2.6.4 from the same book.

    // 1 <= tau1a < tau1b. Controls the alpha jump size during the bracketing phase of
    // the search.
    const double tau1a = 1.4;
    const double tau1b = 9;

    // it must be the case that 0 < tau2 < tau3 <= 1/2 for the algorithm to function
    // correctly but the specific values of tau2 and tau3 aren't super important.
    const double tau2 = 1.0/10.0;
    const double tau3 = 1.0/2.0;


    // Stop right away and return a step size of 0 if the gradient is 0 at the starting point
    if (std::abs(d0) <= std::abs(f0)*std::numeric_limits<double>::epsilon())
        return 0;

    // Stop right away if the current value is good enough according to min_f
    if (f0 <= min_f)
        return 0;

    // Figure out a reasonable upper bound on how large alpha can get.
    const double mu = (min_f - f0) / (rho * d0);

    double alpha = 1;
    if (mu < 0)
        alpha = -alpha;
    alpha = put_in_range(0, 0.65*mu, alpha);


    double last_alpha = 0;
    double last_val = f0;
    double last_val_der = d0;

    // The bracketing stage will find a range of points [a,b]
    // that contains a reasonable solution to the line search
    double a, b;

    // These variables will hold the values and derivatives of f(a) and f(b)
    double a_val, b_val, a_val_der, b_val_der;

    // This thresh value represents the Wolfe curvature condition
    const double thresh = std::abs(sigma*d0);

    unsigned long itr = 0;

    std::vector<double> x_alpha(x.size());

    // do the bracketing stage to find the bracket range [a,b]
    while (true)
    {
        ++itr;
        for (size_t i = 0; i < x.size(); ++i) {
            x_alpha[i] = x[i] + alpha * p[i];
        }
        const double val = f(x_alpha);
        const double val_der = directional_derivative(der(x_alpha),p);

        // we are done with the line search since we found a value smaller
        // than the minimum f value
        if (val <= min_f)
            return alpha;

        if (val > f0 + rho*alpha*d0 || val >= last_val)
        {
            a_val = last_val;
            a_val_der = last_val_der;
            b_val = val;
            b_val_der = val_der;

            a = last_alpha;
            b = alpha;
            break;
        }
        if (std::abs(val_der) <= thresh)
            return alpha;

        // if we are stuck not making progress then quit with the current alpha
        if (last_alpha == alpha || itr >= max_iter)
            return alpha;

        if (val_der >= 0)
        {
            a_val = val;
            a_val_der = val_der;
            b_val = last_val;
            b_val_der = last_val_der;

            a = alpha;
            b = last_alpha;
            break;
        }

        const double temp = alpha;
        // Pick a larger range [first, last].  We will pick the next alpha in that
        // range.
        double first, last;
        if (mu > 0)
        {
            first = std::min(mu, alpha + tau1a*(alpha - last_alpha));
            last  = std::min(mu, alpha + tau1b*(alpha - last_alpha));
        }
        else
        {
            first = std::max(mu, alpha + tau1a*(alpha - last_alpha));
            last  = std::max(mu, alpha + tau1b*(alpha - last_alpha));
        }

        // pick a point between first and last by doing some kind of interpolation
        if (last_alpha < alpha)
            alpha = last_alpha + (alpha-last_alpha)*poly_min_extrap(last_val, last_val_der, val, val_der, 1e10);
        else
            alpha = alpha + (last_alpha-alpha)*poly_min_extrap(val, val_der, last_val, last_val_der, 1e10);

        alpha = put_in_range(first,last,alpha);

        last_alpha = temp;

        last_val = val;
        last_val_der = val_der;

    }//end bracketing while

    // Now do the sectioning phase from 2.6.4
    while (true)
    {
        ++itr;
        double first = a + tau2*(b-a);
        double last = b - tau3*(b-a);
        // use interpolation to pick alpha between first and last
        alpha = a + (b-a)*poly_min_extrap(a_val, a_val_der, b_val, b_val_der);
        alpha = put_in_range(first,last,alpha);

        for (size_t i = 0; i < x.size(); ++i) {
            x_alpha[i] = x[i] + alpha * p[i];
        }
        const double val = f(x_alpha);
        const double val_der = directional_derivative(der(x_alpha),p);

        // we are done with the line search since we found a value smaller
        // than the minimum f value or we ran out of iterations.
        if (val <= min_f || itr >= max_iter)
            return alpha;

        // stop if the interval gets so small that it isn't shrinking any more due to rounding error 
        if (a == first || b == last)
        {
            return b;
        }

        // If alpha has basically become zero then just stop.  Think of it like this,
        // if we take the largest possible alpha step will the objective function
        // change at all?  If not then there isn't any point looking for a better
        // alpha.
        const double max_possible_alpha = std::max(std::abs(a),std::abs(b));
        if (std::abs(max_possible_alpha*d0) <= std::abs(f0)*std::numeric_limits<double>::epsilon())
            return alpha;


        if (val > f0 + rho*alpha*d0 || val >= a_val)
        {
            b = alpha;
            b_val = val;
            b_val_der = val_der;
        }
        else
        {
            if (std::abs(val_der) <= thresh)
                return alpha;

            if ( (b-a)*val_der >= 0)
            {
                b = a;
                b_val = a_val;
                b_val_der = a_val_der;
            }

            a = alpha;
            a_val = val;
            a_val_der = val_der;
        }//end else
    }//end sectioning while
}//end dlib_line_search

// Function to reset the Hessian matrix to the Identity matrix
void reset_hessian(std::vector<std::vector<double>>& H) {
    for (size_t i = 0; i < H.size(); ++i) {
        for (size_t j = 0; j < H[0].size(); ++j) {
            H[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
}

std::vector<double> project_bounds(const std::vector<double>& x, 
                                   const double lower, 
                                   const double upper) {
    if (lower > upper) {
        throw std::runtime_error("Lower bound is greater than upper bound...");
        exit(0);
    }
    std::vector<double> x_projected(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        //x_projected[i] = std::min(std::max(x[i], lower), upper);
        // Ensure that lower bounds are not greater than upper bounds
        if(x[i] < lower) {
            x_projected[i] = lower;
        }
        if (x_projected[i] > upper) {
            x_projected[i] = upper;
        }
    }
    return x_projected;
}

// Limited Memory BFGS-B
std::vector<double> lbfgsb_update(const std::vector<double> &g, std::deque<std::vector<double>> &s_history, 
                                  std::deque<std::vector<double>> &y_history, std::deque<double> &rho_history,
                                  double gamma_k, int m) {
    // Two-loop recursion to compute the search direction 'z'
    std::vector<double> q = g;
    std::vector<double> alphas(m, 0.2);

    // First loop
    for (int i = s_history.size() - 1; i >= 0; --i) {
        alphas[i] = rho_history[i] * dot_product(s_history[i], q);
        q = subtract_vectors(q, scale_vector(y_history[i], alphas[i]));
    }

    std::vector<double> p = scale_vector(q, gamma_k);

    // Second loop
    for (int i = 0; i < s_history.size(); ++i) {
        double beta = rho_history[i] * dot_product(y_history[i], p);
        p = add_vectors(p, scale_vector(s_history[i], alphas[i] - beta));
    }

    return scale_vector(p, -1.0);  // descent direction
}

std::pair<std::vector<double>, std::vector<double>> lbfgsb_step(std::function<double(std::vector<double> &)> func,double alpha, std::vector<double> p,std::string algorithm,
    const std::vector<double>& x, const std::vector<double>& g,
    const double lower, const double upper, // Added bounds as pair
    std::deque<std::vector<double>>& s_history, std::deque<std::vector<double>>& y_history,
    std::deque<double>& rho_history, double& gamma_k, int m) {
    
    // Take a step along the search direction z
    std::vector<double> x_new(x.size());
    for (int j = 0; j < x.size(); j++) {
        x_new[j] = x[j] + alpha * p[j];
    }
    std::vector<double> x_new_unbounded = x;
    if(algorithm == "lbfgsb") {
        try {
            // Your existing code where you call project_bounds
            x_new = project_bounds(x_new_unbounded, lower, upper);
            //std::vector<double> x_projected = project_bounds(x, lower_bounds, upper_bounds);
        } catch (const std::runtime_error& e) {
            std::cerr << "Error in project_bounds: " << e.what() << std::endl;
            // Handle the error, perhaps exit the function or attempt recovery
        }
    }//end boxed LBFGS
    
    // Compute the new gradient
    std::vector<double> g_new = rosenbrock_gradient_regularized(x_new);
    
    // Update L-BFGS-B history
    std::vector<double> s_new = subtract_vectors(x_new, x);
    std::vector<double> y_new = subtract_vectors(g_new, g);
    double rho_new = 1.0 / dot_product(y_new, s_new);
    gamma_k = dot_product(s_new, y_new) / dot_product(y_new, y_new);

    if (s_history.size() == m) {
        s_history.pop_front();
        y_history.pop_front();
        rho_history.pop_front();
    }
    s_history.push_back(s_new);
    y_history.push_back(y_new);
    rho_history.push_back(rho_new);

    return {x_new, g_new};
}

/* The Broyden–Fletcher–Goldfarb–Shanno update */
void bfgs_update(std::vector<std::vector<double>>& H,
                 std::vector<double> delta_x,
                 std::vector<double> delta_g, 
                 double delta_dot) {  

    double yt_s = 1.0 / dot_product(delta_g, delta_x); // y^T s

    if (std::abs(yt_s) < 1e-10) {
        std::cout << "yt_s is too close to zero, resetting H to identity matrix: " << yt_s << std::endl;
        int n = delta_x.size();
        H.clear();
        H.resize(n, std::vector<double>(n, 0));
        for (int i = 0; i < n; i++) {
            H[i][i] = 1.0; // Reset H to identity matrix
        }
        return;
    }

    int n = delta_x.size();
    std::vector<std::vector<double>> I(n, std::vector<double>(n, 0)); // Identity matrix
    for (int i = 0; i < n; i++) {
        I[i][i] = 1;
    }

    std::vector<std::vector<double>> term1 = outer_product(delta_x, delta_g); // s y^T
    std::vector<std::vector<double>> term2 = outer_product(delta_g, delta_x); // y s^T
    std::vector<std::vector<double>> term3 = outer_product(delta_x, delta_x); // s s^T

    for (int row = 0; row < H.size(); row++) {
        for (int col = 0; col < H[row].size(); col++) {
            H[row][col] = (I[row][col] - term1[row][col] * yt_s) * H[row][col] * (I[row][col] - term2[row][col] * yt_s) + term3[row][col] * yt_s;
        }
    }
}


template <typename Vector>
void enhanced_bfgs_update(std::vector<std::vector<double>>& H, 
                          const Vector& x,
                          const Vector& funct_derivative,
                          const Vector& prev_x,
                          const Vector& prev_derivative,
                          bool& been_used,
                          bool& been_used_twice) 
{
    if (!been_used) {
        // Set H to identity matrix
        size_t n = x.size();
        H.clear();
        H.resize(n, std::vector<double>(n, 0));
        for (size_t i = 0; i < n; i++) {
            H[i][i] = 1.0;
        }
        been_used = true;
    } else {
        Vector delta = subtract_vectors(x,prev_x);
        Vector gamma = subtract_vectors(funct_derivative,prev_derivative);

        double dg = dot_product(delta, gamma);

        // Initial H adjustment (as in Code 1)
        if (!been_used_twice) {
            double gg = dot_product(gamma, gamma);
            if (std::abs(gg) > std::numeric_limits<double>::epsilon()) {
                double temp = std::max(0.01, std::min(100.0, dg / gg));
                H.clear();
                H.resize(x.size(), std::vector<double>(x.size(), temp));
                been_used_twice = true;
            }
        }

        double yt_s = 1.0 / dg;
        if (std::abs(yt_s) > 1e-15) {
            // Matrix operations for BFGS update
            for (size_t i = 0; i < x.size(); i++) {
                for (size_t j = 0; j < x.size(); j++) {
                    double syT = delta[i] * gamma[j];
                    double ysT = gamma[i] * delta[j];
                    double ssT = delta[i] * delta[j];
                    H[i][j] = (1 - syT * yt_s) * H[i][j] * (1 - ysT * yt_s) + ssT * yt_s;
                }
            }
        } else {
            // Fallback to identity matrix
            size_t n = x.size();
            H.clear();
            H.resize(n, std::vector<double>(n, 0));
            for (size_t i = 0; i < n; i++) {
                H[i][i] = 1.0;
            }
            been_used_twice = false;
        }
    }
}



void dfp_update(std::vector<std::vector<double>>& H, std::vector<double> delta_x, std::vector<double> delta_g) {
    double term2 = dot_product(delta_x, delta_g); // O(n)
    if (std::abs(term2) < 1e-10) {
        //std::cout << "term2 is too close to zero: " << term2 << std::endl;
        return;
    }

    std::vector<std::vector<double>> term1 = outer_product(delta_x, delta_x); // O(n^2)

    // Simplified term: denom = H_k^T delta_g delta_g^T
    // Calculate H * delta_g
    std::vector<double> H_delta_g = matvec_product(H, delta_g); // O(n^2)

    // Calculate delta_g^T * H * delta_g
    double term5 = dot_product(delta_g, H_delta_g); // O(n)
    if (std::abs(term5) < 1e-15) {
        std::cout << "term5 is too close to zero: " << term5 << std::endl;
        return;
    }

    for (int row = 0; row < H.size(); row++) {
        for (int col = 0; col <= row; col++) { // Symmetry saves about half the calculations
            double denom_term = H_delta_g[row] * delta_g[col] + H_delta_g[col] * delta_g[row]; // 2 * (H * delta_g)_row * delta_g_col
            H[row][col] = H[row][col] + term1[row][col] / term2 - denom_term / term5;
            if (row != col) {
                H[col][row] = H[row][col];  // H is symmetric
            }
        }
    }
}

double optimize(std::function<double(std::vector<double> &)> func, std::function<std::vector<double>(const std::vector<double>&)> der, std::vector<double> x0,std::string algorithm,const double tol,const int max_iter,const double lower, const double upper) {
    double min_value = std::numeric_limits<double>::max();
    double alpha = 1.0,max_alpha = 1.0, lambda = 15.0;
    double min_alpha, fnew, delta_dot;
    double current_fun, previous_fun, min_delta = 1e-25;
    
    std::vector<double> x = x0;
    std::vector<double> p, g, delta_x, delta_g, x_new;
    std::vector<double> x_new_a(x.size());
    if(algorithm.find("dfp") != std::string::npos) {
        min_alpha = 1e-6;
    } else {
        min_alpha = 1e-4;
    }

    std::pair<std::vector<double>, std::vector<double>> result;
    // L-BFGS-B init
    std::deque<std::vector<double>> s_history, y_history;
    std::deque<double> rho_history;
    double gamma_k = 1.0;  // initial scaling factor
    int m = 3;  // history size

    // Initialize the Hessian matrix to identity matrix
    std::vector<std::vector<double>> H(x0.size(), std::vector<double>(x0.size(), 0));
    for (int i = 0; i < x0.size(); i++) {
        H[i][i] = 1; // along the diagonals, place 1's to create Identity Matrix
    }//end for
    
    int early_stopping = 0;
    int i;

    bool been_used = false;
    bool been_used_twice = false;

    // Main loop
    for (i = 0; i < max_iter; i++) {
        //current_fun = func(x);
        g = der(x); // Compute Gradient
        /*for(int j=0;j<g.size();j++){
            if(!isValidDouble(g[j])) {
                std::cout << "i: "<<i << "j: "<<j<<"\nnonvalid value, exiting..."<<std::endl;
                exit(0);
            }
        }
        */
        // Check if the length of gradient vector is less than our tolerance
        if (norm(g) < tol) {
            min_value = std::min(min_value, func(x));
            global_min = min_value;
            if(global_min > 0.0 && min_value > 0.0) { continue;}
            std::cout << i <<"# it" << std::endl;
            std::cout << "\nnorm(g): New Global Minimum: " << global_min << " with parameters:" <<std::endl;
            best_params = {};
            for (int i=0;i<=x.size();i++){
                best_params.push_back(x[i]);
                std::cout<< "x["<<i<<"]: " << best_params[i]<<std::endl;
            }
            return min_value;
        }//end if

        // Compute Search Direction
        p = matvec_product(H, g);
        p = scale_vector(p, -1.0); //opposite of greatest increase

        /*** Calculate optimal step size in the search direction p ***/
        alpha = dlib_line_search(func,der, func(x),dot_product(g, p), // Current derivative value along p
            0.15,                  // rho,  BFGS default=0.01, W&PQ=0.15, R=0.25, ALL=0.015
            0.90,                  // sigma, BFGS default=0.90, W&PQ=0.90, R=0.75, ALL=0.90
            1e-9,                // min_f, 1e-12 to 1e-7 to 1e-9
            100,x,p,lambda
        );

        if(!isValidDouble(alpha)) {
            std::cout << "step size can't be nan..." << std::endl;
            alpha = 1.0;
        }
        if (alpha < min_alpha) {
            alpha = min_alpha;
        } else if (alpha > max_alpha) {
            alpha = max_alpha;
        }

        if(algorithm.find("lbfgs") != std::string::npos) {  
            result = lbfgsb_step(func,alpha, p,algorithm,x0,g,lower,upper,s_history, y_history, rho_history, gamma_k, m);
            x = result.first;
            g = result.second;
        } else {
            x_new = add_vectors(x, scale_vector(p, alpha)); // Update the current point x by taking a step of size alpha in the direction p.
            delta_x = subtract_vectors(x_new, x);           // Compute the difference between the new point and the old point, delta_x
            delta_g = subtract_vectors(der(x_new), g);      // Compute the difference in the gradient at the new point and the old point, delta_g.
            delta_dot = dot_product(delta_g, delta_x);
            if (std::abs(delta_dot) < 1e-10) {
                //printVectors(i, fnew, x, delta_x, delta_g, alpha);
                //break;
                if(alpha == min_alpha){
                    alpha *= 2.5;
                } else {
                    alpha = 0.5;
                }
                //std::cout << "new alpha = " << alpha << std::endl;
                std::vector<double> new_dir = scale_vector(p, alpha);
                add_vectors(x, new_dir);
                if(early_stopping < 10) {
                    early_stopping++;
                } else {
                    break;
                }
            }
            if (algorithm == "bfgs") { // Update the inverse Hessian approximation using BFGS             
                bfgs_update(H, delta_x, delta_g, delta_dot);
                // INSERT IT HERREEEE !!!! enhanced_bfgs_update()
                //enhanced_bfgs_update(H, x_new, der(x), x, g, been_used, been_used_twice);
            } else {// Update the approximation of the inverse Hessian using DFP
                dfp_update(H, delta_x, delta_g);
            }
            /*if (std::abs(current_fun - previous_fun) < min_delta)
            {
                std::cout << "delta fun < " << min_delta;
                break;
            }
            previous_fun = current_fun;
            */
            x = x_new;
            min_value = std::min(min_value, func(x));

        }//end else
        if(max_iter == i) {
            std::cout << "Max iterations reached"<<std::endl;
        }
    }// end main loop
    
    std::cout << i << " iterations" << std::endl;
    print_vector(x);
    min_value = std::min(min_value, func(x));
    return min_value; 
}// end optimize

double minimize(std::function<double(std::vector<double> &)> func,std::function<std::vector<double>(const std::vector<double>&)> der,  std::vector<double> x0, std::string name, 
             const int pop_size,const int max_gens,const int dim, std::string algorithm, const double lower, const double upper) {
    //global_min = std::numeric_limits<double>::max();
    double minima = optimize(func,der, x0, algorithm, 1e-15, 2000, lower, upper);
    std::cout << "Predicted Global minimum for " << name << ": " << minima <<std::endl;
    return minima;
}
