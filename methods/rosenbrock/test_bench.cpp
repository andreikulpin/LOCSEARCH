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
    searchMethod.getOptions().mHInit = initH;
    searchMethod.getOptions().mDoTracing = false;
    searchMethod.getOptions().mInc = 5.0;
    searchMethod.getOptions().mDec = - 0.5;
    searchMethod.getOptions().maxUnsuccessStepsNumber = 50;
    searchMethod.getOptions().maxStepsNumber = 100000;
    searchMethod.getOptions().mEps = eps;
    
    double x[dim];
    /*for (int i = 0; i < dim; i++) {
        double a = bm->getBounds()[i].first;
        double b = bm->getBounds()[i].second;
        x[i] = (b + a) / 2.0;
    }*/
    snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
    
    std::cout << "*************Testing benchmark**********" << std::endl;
    std::cout << bm->getDesc() << std::endl;
    
    double v;
    searchMethod.search(x, v);
    std::cout << "Found v = " << v << "\n";
    std::cout << "the difference is " << v - bm->getGlobMinY() << std::endl;
    std::cout << " at " << snowgoose::VecUtils::vecPrint(dim, x) << "\n";
    
    int count;
    //std::shared_ptr<COMPI::Functor<double>> objPtr = mpp->mObjectives.at(0).get();
    if (auto obj = dynamic_cast<COMPI::FuncCnt<double> *>(mpp->mObjectives.at(0).get())) {
        count = obj->mCounters.mFuncCalls;
    }
    std::cout << "Number of objective calls is " << count << std::endl;
    
    std::cout << "****************************************" << std::endl << std::endl;
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

