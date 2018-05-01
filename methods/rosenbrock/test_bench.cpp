/* 
 * File:   test_bench.cpp
 * Author: kulpin
 *
 * Created on April 17, 2018, 12:42 AM
 */

#include <iostream>
#include "testfuncs/manydim/benchmarks.hpp"
#include "rosenbrock_method.hpp"
#include <oneobj/contboxconstr/benchmarkfunc.hpp>

using BM = Benchmark<double>;

bool testBench(std::shared_ptr<BM> bm, double eps) {
    const int dim = bm->getDim();
    
    double initH[dim];
    snowgoose::VecUtils::vecSet(dim, .1, initH);
    
    OPTITEST::BenchmarkProblemFactory problemFactory(bm);
    COMPI::MPProblem<double> *mpp = problemFactory.getProblem();
    
    LOCSEARCH::RosenbrockMethod<double> searchMethod(*mpp);
    searchMethod.getOptions().mHInit = initH;
    searchMethod.getOptions().mDoTracing = false;
    searchMethod.getOptions().mInc = 5.0;
    searchMethod.getOptions().mDec = - 0.5;
    searchMethod.getOptions().maxUnsuccessStepsNumber = 50;
    searchMethod.getOptions().maxStepsNumber = 100000;
    searchMethod.getOptions().mEps = eps;
    
    double x[dim];
    
    for (int i = 0; i < dim; i++) {
        double a = bm->getBounds()[i].first;
        double b = bm->getBounds()[i].second;
        x[i] = (b + a) / 2.0;
    }
    
    std::cout << "*************Testing benchmark**********" << std::endl;
    //std::cout << bm;
    
    double v;
    searchMethod.search(x, v);
    std::cout << "v = " << v << std::endl;
    std::cout << "the difference is " << v - bm->getGlobMinY() << std::endl;
    std::cout << "****************************************" << std::endl << std::endl;
}

int main(int argc, char** argv) {
    const int dim = argc > 1 ? atoi(argv[1]) : 50;
    double eps = argc > 2 ? atof(argv[2]) : 0.01;
    
    Benchmarks<double> tests;
    for (auto bm : tests) {
        testBench(bm, eps);
    }
}

