/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   test_acd_bench.cpp
 * Author: kulpin
 *
 * Created on April 27, 2018, 12:48 AM
 */

#include <iostream>
#include <funccnt.hpp>
#include <methods/lins/goldsec/goldsec.hpp>
#include "acd_bench.hpp"
#include <oneobj/contboxconstr/benchmarkfunc.hpp>

using BM = Benchmark<double>;

bool testBench(std::shared_ptr<BM> bm, double eps) {
    const int dim = bm->getDim();
    
    OPTITEST::BenchmarkProblemFactory problemFactory(bm);
    COMPI::MPProblem<double> *mpp = problemFactory.getProblem();
    
    LOCSEARCH::GoldenSecLS<double>* locs = new LOCSEARCH::GoldenSecLS<double>(*mpp);
    locs->getOptions().mSInit = 0.1;
    locs->getOptions().mDelta = 0.02;
    locs->getOptions().mMaxBackSteps = 16;
    locs->getOptions().mDoTracing = false;
    
    LOCSEARCH::AdvancedCoordinateDescentBench<double> searchMethod(*mpp, bm->getGlobMinY());
    searchMethod.getLineSearch().reset(locs);    
    searchMethod.getOptions().mHInit = .1;
    searchMethod.getOptions().mDoTracing = false;
    searchMethod.getOptions().mGradLB = 0;
    searchMethod.getOptions().mSearchType = LOCSEARCH::AdvancedCoordinateDescentBench<double>::SearchTypes::NO_DESCENT;
    searchMethod.getOptions().mVicinityAdaptation = LOCSEARCH::AdvancedCoordinateDescentBench<double>::UNIFORM_ADAPTATION;
    searchMethod.getOptions().maxStepsNumber = 100000;
    searchMethod.getOptions().mEps = eps;
    
    double x[dim];
    snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
    
    double v;
    
    /*
    std::cout << "*************Testing benchmark**********" << std::endl;
    std::cout << bm->getDesc() << std::endl;
    
    searchMethod.search(x, v);
    std::cout << "v = " << v << std::endl;
    std::cout << "the difference is " << v - bm->getGlobMinY() << std::endl;
    std::cout << " at " << snowgoose::VecUtils::vecPrint(dim, x) << "\n";
    
    int count;
    if (auto obj = dynamic_cast<COMPI::FuncCnt<double> *>(mpp->mObjectives.at(0).get())) {
        count = obj->mCounters.mFuncCalls;
    }
    std::cout << "Number of objective calls is " << count << std::endl;
     * */
    
    
    searchMethod.search(x, v);
    
    if (SGABS(v - bm->getGlobMinY()) < 0.1) {
        std::cout << bm->getDesc() << "\t";
        std::cout << bm->getGlobMinY() << "\t";
        std::cout << v << "\t";

        int count;
        if (auto obj = dynamic_cast<COMPI::FuncCnt<double> *>(mpp->mObjectives.at(0).get())) {
            count = obj->mCounters.mFuncCalls;
        }
        std::cout << count << std::endl;
    }
}

int main(int argc, char** argv) {
    const int dim = argc > 1 ? atoi(argv[1]) : 50;
    double eps = argc > 2 ? atof(argv[2]) : 0.01;
    
    /*auto bm = std::make_shared<RosenbrockBenchmark<double>>(dim);
    testBench(bm, eps);*/
    
    Benchmarks<double> tests;
    for (auto bm : tests) {
        testBench(bm, eps);
    }
    return 0;
}

