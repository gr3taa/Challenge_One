#ifndef MET_HPP
#define MET_HPP

#include"helper.hpp"
#include<iostream>
#include<cmath>

// this class is used to set the method to compute alpha
enum class ALPHA_STRATEGY{
    Armijo,
    Inverse_decay,
    Exponential_decay
};

// this class is used to set the method to compute the gradient
enum class GRAD_METHOD{
    Exact_gradient,
    Finite_differences
};

//definition of the exact gradient of function to minimize
vector<double> grad(vector<double> c){
    vector<double> ret;
    ret.resize(2);
    ret[0]=(c[1] + 16*pow(c[0],3) + 3);
    ret[1]=( c[0] + 2 * c[1]);
    return ret;
}

//definition of the function to minimize
double func(vector<double> c){
    return 3*c[0] + c[0]*c[1] + 4*pow(c[0], 4) + pow(c[1],2);
}


struct parameters
{
    function<double(vector<double>)> f = func;
    function<vector<double>(vector<double>)> grad_f = grad;
    vector<double> init_cond = {0 , 0};//initial condition
    double tol_s ;
    double tol_r;
    double alpha_zero ;
    int max_it;
    double sigma ; // parameter for Armijo rule
    double mu ; 
    double eta ;
    ALPHA_STRATEGY strategy;
    GRAD_METHOD g;
};

//this function prints the parameters
void print_parameters(parameters p) {
    cout<<"Parameters used:"<<endl;
    cout<<"The Initial Condition x0 = [ " << p.init_cond[0] <<", "
         << p.init_cond[1] << "]"<<endl;
    cout<<"Eps_s = "<<p.tol_s<<endl;
    cout<<"Eps_r = "<<p.tol_r<<endl;
    cout<<"Alpha_0 = "<<p.alpha_zero<<endl;
    cout<<"Max iteration = "<<p.max_it<<endl;
    cout<<"Sigma = "<<p.sigma<<endl;
    cout<<"Mu = "<<p.mu<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------"<<endl;
}

//this function prints the methods used
void print_methods(parameters p, const vector<double> &x_new, int it){
    cout<<"Strategy used for computation of alpha_k is: ";
    if(p.strategy == ALPHA_STRATEGY::Armijo)
        cout<<"Armijo rule"<<endl;
    else if(p.strategy == ALPHA_STRATEGY::Inverse_decay)
        cout<<"Inverse decay"<<endl;
    else if(p.strategy == ALPHA_STRATEGY::Exponential_decay)
        cout<<"Exponential decay"<<endl;
    cout<<"Method used for the computation of the gradient: ";
    if(p.g == GRAD_METHOD::Exact_gradient)
        cout<<"Exact gradient"<<endl;
    else if(p.g == GRAD_METHOD::Finite_differences)
        cout<<"Finite differences"<<endl;
    cout<<"The method converges in "<< it<<" iterations"<<endl;
    cout<<"The minimum is x = [ "<<x_new[0]<<", "<< x_new[1]<<"]\n"<<endl;
    cout<<"------------------------------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"\n";
}

// this function compute alpha through the Armijo rule
double Armijo_rule(parameters p, const vector<double>& x, double alpha_init, double par) {
    while(p.f(x)-p.f(x-alpha_init*p.grad_f(x)) < par * alpha_init * vec_norm(p.grad_f(x)) * vec_norm(p.grad_f(x)))           
        alpha_init = alpha_init/2;
    return alpha_init;
}


// this function compute alpha through the exponential decay
double Exponential(double alpha_init, double mu, int it) {
    return alpha_init * exp(-(mu * it));
}


// this function compute alpha through the inverse decay
double Inverse(double alpha_init, double mu, int it) {
    return alpha_init/(1 + mu * it);
}

//through a template this function choose what method for the computation of the gradient has been required
template<GRAD_METHOD G>
vector<double> compute_gradient(vector<double> v, parameters p){
    if constexpr (G == GRAD_METHOD::Exact_gradient)
        return p.grad_f(v);
    if constexpr (G == GRAD_METHOD::Finite_differences)
        return diff_finite(p.f, 0.0001, v);
}

//this functions solve the minimization problem through the gradient method
template<ALPHA_STRATEGY C, GRAD_METHOD G>
vector<double> grad_method( parameters p){
    p.strategy = C;
    p.g = G;
    vector<double> x_vec = p.init_cond; // x_0
    vector<double> x_new = x_vec -  p.alpha_zero * compute_gradient<G>(x_vec, p); // x_1
    double alpha = p.alpha_zero;
    int it = 1;
    //check the convergence conditions
    while(vec_norm(x_new - x_vec) > p.tol_s && vec_norm(compute_gradient<G>(x_vec, p)) > p.tol_r && it < p.max_it){
        ++it;
        x_vec = x_new; //x_k
        //choose the method to compute alpha
        if constexpr (C == ALPHA_STRATEGY::Armijo){
            alpha = Armijo_rule(p, x_vec, p.alpha_zero, p.sigma);
        } else if constexpr (C == ALPHA_STRATEGY::Exponential_decay){
            alpha = Exponential(p.alpha_zero, p.mu, it);
        } else if constexpr (C == ALPHA_STRATEGY::Inverse_decay){
            alpha = Inverse(p.alpha_zero, p.mu, it);
        } else {
            throw invalid_argument("Invalid alpha strategy");
        }
        x_new = x_vec -  alpha * compute_gradient<G>(x_vec, p); //x_(k+1)
    }
    cout<<"The gradient method\n";
    print_methods(p, x_new, it);
    return x_new;
}


//this functions solve the minimization problem through the heavy ball method
template<ALPHA_STRATEGY C, GRAD_METHOD G>
vector<double> heavy_ball_method(parameters p){
    p.strategy = C;
    p.g = G;
    vector<double> x_vec = p.init_cond; //x_0
    vector<double> x = p.init_cond - p.alpha_zero * compute_gradient<G>(p.init_cond, p); //x_1
    vector<double> x_new = x;
    double alpha = p.alpha_zero;
    int it = 1;
    //check the convergence conditions
    while(vec_norm(x - x_vec) > p.tol_s && vec_norm(compute_gradient<G>(x_vec, p)) > p.tol_r && it < p.max_it){
        ++it;
        //choose the method to compute alpha
        if constexpr (C == ALPHA_STRATEGY::Armijo){
            cerr<<"It is not possible to use Armijo rule to compute alpha\n";
            exit(1);
        } else if constexpr (C == ALPHA_STRATEGY::Exponential_decay){
            alpha = Exponential(p.alpha_zero, p.mu, it);
        } else if constexpr (C == ALPHA_STRATEGY::Inverse_decay){
            alpha = Inverse(p.alpha_zero, p.mu, it);
        } else {
            throw invalid_argument("Invalid alpha strategy");
        }
        x_new = x - alpha * compute_gradient<G>(x, p) + p.eta * (x - x_vec); //x_(k+1)
        x_vec = x; //x_(k-1)
        x = x_new; //x_k
    }
    cout<<"The heavy ball method\n";
    print_methods(p, x_new, it);
    return x_new;
}


//this functions solve the minimization problem through the Nesterov method
template<ALPHA_STRATEGY C, GRAD_METHOD G>
vector<double> Nesterov_method(parameters p){
    p.strategy = C;
    p.g = G;
    vector<double> x_vec = p.init_cond; //x_0
    vector<double> x = p.init_cond - p.alpha_zero * compute_gradient<G>(p.init_cond, p); //x_1
    vector<double> x_new = x;
    x_new.resize(x.size());
    double alpha = p.alpha_zero;
    int it = 1;
    //check the convergence conditions
    while(vec_norm(x - x_vec) > p.tol_s && vec_norm(compute_gradient<G>(x_vec, p)) > p.tol_r && it < p.max_it){
        ++it;
        //choose the method to compute alpha
        if constexpr (C == ALPHA_STRATEGY::Armijo){
            cerr<<"It is not possible to use Armijo rule to compute alpha\n";
            exit(1);
        } else if constexpr (C == ALPHA_STRATEGY::Exponential_decay){
            alpha = Exponential(p.alpha_zero, p.mu, it);
        } else if constexpr (C == ALPHA_STRATEGY::Inverse_decay){
            alpha = Inverse(p.alpha_zero, p.mu, it);
        } else {
            throw invalid_argument("Invalid alpha strategy");
        }
        vector<double> y = x + p.eta * (x - x_vec);
        x_new = y - alpha * compute_gradient<G>(y, p) ;//x_(k+1)
        x_vec = x;//x_(k-1)
        x = x_new;//x_k
    }
    cout<<"The Nesterov method\n";
    print_methods(p, x_new, it);
    return x_new;
}












#endif