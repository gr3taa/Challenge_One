#include<vector>
#include<cmath>
#include<functional>
#include<iostream>
#include"met.hpp"
#include"GetPot"

using namespace std;

int main(int argc, char ** argv){
    parameters p; 

    GetPot c_line(argc, argv);
    p.mu = c_line("mu", 0.3);
    p.alpha_zero = c_line("alpha_zero", 0.1);
    p.eta = c_line("eta", 0.9);
    p.max_it = c_line("max_it", 1000);
    p.sigma = c_line("sigma", 0.3);
    p.tol_r = c_line("tol_r", 1e-6);
    p.tol_s = c_line("tol_s", 1e-6);

    print_parameters(p);

    vector<double> grad_res = grad_method<ALPHA_STRATEGY::Armijo, GRAD_METHOD::Exact_gradient>(p);  
    vector<double> grad_h_b = heavy_ball_method<ALPHA_STRATEGY::Inverse_decay, GRAD_METHOD::Exact_gradient>(p);
    vector<double> grad_Nesterov = Nesterov_method<ALPHA_STRATEGY::Inverse_decay, GRAD_METHOD::Finite_differences>(p);
    
    return 0;
}





