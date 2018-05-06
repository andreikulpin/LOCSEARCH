/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include <iostream>
#include <funccnt.hpp>
#include "testfuncs/manydim/benchmarks.hpp"
#include <oneobj/contboxconstr/benchmarkfunc.hpp>
#include "../rosenbrock/rosenbrock_bench.hpp"
#include "../rosenbrock/rosenbrock_bench_old.hpp"
#include <cmath>

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
    const int dim = 10;
    
    int minPow = 4;
    int maxPow = 16;
    int stepPow = 2;
    
    for (int i = minPow; i <= maxPow; i += stepPow) {
        const double eps = std::pow(10, -i);

        auto bm = std::make_shared<RosenbrockBenchmark<double>>(dim);
        
        std::cout << eps << "\t";
        testRosenbrock(bm, eps);
        std::cout << "\t";
        testRosenbrockOld(bm, eps);
        std::cout << "\n";
    }
}
