/* 
 * File:   test_bench.cpp
 * Author: kulpin
 *
 * Created on April 17, 2018, 12:42 AM
 */

#include <iostream>
#include <funccnt.hpp>
#include "testfuncs/manydim/benchmarks.hpp"
#include "rosenbrock_bench.hpp"
#include <oneobj/contboxconstr/benchmarkfunc.hpp>

using BM = Benchmark<double>;

bool testBench(std::shared_ptr<BM> bm, double eps) {
    const int dim = bm->getDim();
    
    double initH[dim];
    snowgoose::VecUtils::vecSet(dim, .1, initH);
    
    OPTITEST::BenchmarkProblemFactory problemFactory(bm);
    COMPI::MPProblem<double> *mpp = problemFactory.getProblem();
    //auto obj = dynamic_cast<std::shared_ptr<COMPI::FuncCnt<double>>>(mpp->mObjectives.at(0));
    //auto obj = std::make_shared<COMPI::FuncCnt<double>>(objPtr);
    
    LOCSEARCH::RosenbrockBenchMethod<double> searchMethod(*mpp, bm->getGlobMinY());
    searchMethod.getOptions().mHInit = std::vector<double>(initH, initH + dim);
    searchMethod.getOptions().mDoTracing = false;
    searchMethod.getOptions().mInc = 5.0;
    searchMethod.getOptions().mDec = 0.5;
    searchMethod.getOptions().mMaxStepsNumber = 100000;
    searchMethod.getOptions().mHLB = searchMethod.getOptions().mMinGrad * 1e-2;
    searchMethod.getOptions().mEps = eps;
    
    double x[dim];
    /*for (int i = 0; i < dim; i++) {
        double a = bm->getBounds()[i].first;
        double b = bm->getBounds()[i].second;
        x[i] = (b + a) / 2.0;
    }*/
    snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
    
    double v;
    
    searchMethod.search(x, v);
    
    std::cout << bm->getDesc() << "\t";
    std::cout << bm->getGlobMinY() << "\t";
    std::cout << v << "\t";

    int count;
    if (auto obj = dynamic_cast<COMPI::FuncCnt<double> *>(mpp->mObjectives.at(0).get())) {
        count = obj->mCounters.mFuncCalls;
    }
    std::cout << count << std::endl;
}

int main(int argc, char** argv) {
    const int dim = argc > 1 ? atoi(argv[1]) : 50;
    const double eps = argc > 2 ? atof(argv[2]) : 0.01;
    
    /*auto bm = std::make_shared<RosenbrockBenchmark<double>>(dim);
    testBench(bm, eps);*/
    
    Benchmarks<double> tests;
    for (auto bm : tests) {
        testBench(bm, eps);
    }
}

