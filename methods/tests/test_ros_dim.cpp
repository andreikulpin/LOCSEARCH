/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   test_ros_dim.cpp
 * Author: andrei
 *
 * Created on May 7, 2018, 12:34 AM
 */

#include <iostream>
#include <funccnt.hpp>
#include "testfuncs/manydim/benchmarks.hpp"
#include <oneobj/contboxconstr/benchmarkfunc.hpp>
#include "../rosenbrock/rosenbrock_bench.hpp"
#include "../rosenbrock/rosenbrock_bench_old.hpp"

using BM = Benchmark<double>;

bool testRosenbrock(std::shared_ptr<BM> bm, double eps) {
    const int dim = bm->getDim();
    
    double initH[dim];
    snowgoose::VecUtils::vecSet(dim, .1, initH);
    
    OPTITEST::BenchmarkProblemFactory problemFactory(bm);
    COMPI::MPProblem<double> *mpp = problemFactory.getProblem();
    
    LOCSEARCH::RosenbrockBenchMethod<double> searchMethod(*mpp, bm->getGlobMinY());
    searchMethod.getOptions().mHInit = std::vector<double>(initH, initH + dim);
    searchMethod.getOptions().mDoTracing = false;
    searchMethod.getOptions().mInc = 5.0;
    searchMethod.getOptions().mDec = 0.5;
    searchMethod.getOptions().mMaxStepsNumber = 500000;
    searchMethod.getOptions().mMinGrad = eps;
    searchMethod.getOptions().mHLB = eps;
    searchMethod.getOptions().mEps = eps;
     
    double x[dim];
    snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
    
    double v;
    
    searchMethod.search(x, v);
    
    auto obj = dynamic_cast<COMPI::FuncCnt<double> *>(mpp->mObjectives.at(0).get());
    int count = obj->mCounters.mFuncCalls;
    std::cout << count;
}

bool testRosenbrockOld(std::shared_ptr<BM> bm, double eps) {
    const int dim = bm->getDim();
    
    double initH[dim];
    snowgoose::VecUtils::vecSet(dim, .1, initH);
    
    OPTITEST::BenchmarkProblemFactory problemFactory(bm);
    COMPI::MPProblem<double> *mpp = problemFactory.getProblem();
    
    LOCSEARCH::RosenbrockBenchMethodOld<double> searchMethod(*mpp, bm->getGlobMinY());
    searchMethod.getOptions().mHInit = initH;
    searchMethod.getOptions().mDoTracing = false;
    searchMethod.getOptions().mInc = 5.0;
    searchMethod.getOptions().mDec = 0.5;
    searchMethod.getOptions().maxUnsuccessStepsNumber = 50;
    searchMethod.getOptions().maxStepsNumber = 500000;
    searchMethod.getOptions().mEps = eps;
    
    double x[dim];
    snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
    
    double v;
    
    bool rv = searchMethod.search(x, v);
    
    auto obj = dynamic_cast<COMPI::FuncCnt<double> *>(mpp->mObjectives.at(0).get());
    int count = obj->mCounters.mFuncCalls;
    std::cout << count;
}

int main(int argc, char** argv) {
    //const int dim = argc > 1 ? atoi(argv[1]) : 50;
    const double eps = 1e-8;
    
    int minN = 1;
    int maxN = 20;
    int stepN = 1;
    
    for (int i = minN; i <= maxN; i += stepN) {
        const int dim = i;
        
        auto bm = std::make_shared<RosenbrockBenchmark<double>>(dim);
        
        std::cout << dim << "\t";
        testRosenbrock(bm, eps);
        std::cout << "\t";
        testRosenbrockOld(bm, eps);
        std::cout << "\n";
    }
}
