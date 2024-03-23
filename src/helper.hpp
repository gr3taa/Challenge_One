#ifndef HELPER_HPP
#define HELPER_HPP


#include<vector>
#include<cmath>
#include<functional>
#include<iostream>
using namespace std;


// this function computes the norm of a standard vector
// it takes in input a vector and returns its norm
double vec_norm(const vector<double>& v){
    size_t n = v.size(); 
    double ret = 0.;
    for(size_t i = 0; i < n; ++i)
        ret += v[i] * v[i];
    return sqrt(ret);
}


//this function computes the product between a scalar and a vector
vector<double> operator*(const double & c, const vector<double>& v){
    vector<double> ret;
    ret.resize(v.size());
    for(size_t i = 0; i < v.size(); ++i)
        ret[i] = c * v[i];
    return ret;
}


//this function computes the difference between two vectors
vector<double> operator-(const vector<double> & lhs,const vector<double> & rhs){
    vector<double> ret;
    if( lhs.size() != rhs.size()){
        cerr<<"Dimensions of the two vectors are different"<<endl;
        return ret;
    }
    ret.resize(lhs.size());
    for(size_t i = 0; i < lhs.size(); ++i)
        ret[i] = lhs[i] - rhs[i];
    return ret;
}


//this function computes the sum between two vectors
vector<double> operator+(const vector<double> & lhs,const vector<double> & rhs){
    vector<double> ret;
    if( lhs.size() != rhs.size()){
        cerr<<"Dimensions of the two vectors are different"<<endl;
        return ret;
    }
    ret.resize(lhs.size());
    for(size_t i = 0; i < lhs.size(); ++i)
        ret[i] = lhs[i] + rhs[i];
    return ret;
}


//this function computes the gradient by finite differences
vector<double> diff_finite (function<double(vector<double>)> f, double h, const vector<double>& v){
    vector<double> ret;
    ret.resize(v.size());
    for (size_t i = 0; i < v.size(); ++i){
        vector<double> v_p = v;
        vector<double> v_m = v;
        v_p[i] += h;
        v_m[i] -= h;
        ret[i] = (f(v_p)-f(v_m))/(2.*h); 
    }
    return ret;
}




#endif