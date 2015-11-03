/* 
 * File:   casadi.hpp
 * Author: Abuenameh
 *
 * Created on 06 November 2014, 17:45
 */

#ifndef CASADI_HPP
#define	CASADI_HPP

#include <typeinfo>

#include <casadi/casadi.hpp>
#include <casadi/solvers/rk_integrator.hpp>
#include <casadi/solvers/collocation_integrator.hpp>
#include <casadi/interfaces/sundials/cvodes_interface.hpp>
#include <casadi/core/function/custom_function.hpp>

using namespace casadi;

#include <boost/date_time.hpp>

using namespace boost::posix_time;

#include <nlopt.hpp>

using namespace nlopt;

#include "gutzwiller.hpp"

struct results {
    double tau;
    double Ei;
    double Ef;
    double Q;
    double p;
    //    vector<vector<double>> bs;
    double U0;
    vector<double> J0;
    vector<complex<double>> b0;
    vector<complex<double>> bf;
    vector<vector<complex<double>>> f0;
    vector<vector<complex<double>>> ff;
    string runtime;
};

class DynamicsProblem {
public:
    DynamicsProblem(double Wi, double Wf, double mu, vector<double> xi);

    static void setup(double W, double mu, vector<double>& xi);

    void evolve(double tau, results& res);

private:

    static SX energy(SX& fin, SX& J, SX& U0, SX& dU, double mu);
    static SX energy(int i, int n, SX& fin, SX& J, SX& U0, SX& dU, double mu);
    static SX canonical(SX& fin, SX& J, SX& U0, SX& dU, double mu);
    static SX canonical(int i, int n, SX& f, SX& J, SX& U0, SX& dU, double mu);

    static double U00;
    static vector<double> J0;

    SXFunction E0;
    
    SXFunction ode_func;
    Integrator integrator;

    static vector<double> x0;

    static vector<vector<complex<double>>> f0;
    vector<vector<complex<double>>> ff;

};

#endif	/* CASADI_HPP */

