/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   test_unimodal.cpp
 * Author: andrei
 *
 * Created on May 5, 2018, 3:38 PM
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>

#include <funccnt.hpp>
#include "testfuncs/manydim/benchmarks.hpp"
#include <oneobj/contboxconstr/benchmarkfunc.hpp>
#include "../advancedcoordescent/acd_bench.hpp"
#include "../lins/goldsec/goldsec.hpp"
#include "../rosenbrock/rosenbrock_bench.hpp"

using BM = Benchmark<double>;

void testBench(std::shared_ptr<BM> bm, double eps, std::ofstream& file);
bool testAcdBench(std::shared_ptr<BM> bm, double eps, std::ofstream& file);
bool testRosenbrockBench(std::shared_ptr<BM> bm, double eps, std::ofstream& file);

int main(int argc, char** argv) {
    std::ofstream file;
    file.open("test_unimodal.out");
    
    const int dim = argc > 1 ? atoi(argv[1]) : 50;
    const double eps = argc > 2 ? atof(argv[2]) : 0.01;
    
    /*auto bm = std::make_shared<RosenbrockBenchmark<double>>(dim);
    testBench(bm, eps, file);*/
    
    Benchmarks<double> tests;
    for (auto bm : tests) {
        if (bm->getMulMod()) {
            testBench(bm, eps, file);
        }
    }
    
    file.close();
    return 0;
}

void testBench(std::shared_ptr<BM> bm, double eps, std::ofstream& file) {
    const int dim = bm->getDim();
    
    file << std::boolalpha; 
    file << bm->getDesc() << "\t";
    file << bm->getGlobMinY() << "\t";
    
    testAcdBench(bm, eps, file);
    file << "\t";
    testRosenbrockBench(bm, eps, file);
    
    
    file << std::endl;
}

bool testAcdBench(std::shared_ptr<BM> bm, double eps, std::ofstream& file) {
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
    
    searchMethod.search(x, v);
    
    file << v << "\t";
    file << std::abs(v - bm->getGlobMinY()) << "\t";
    
    bool isClose = std::abs(v - bm->getGlobMinY()) < eps;
    file << isClose << "\t";

    auto obj = dynamic_cast<COMPI::FuncCnt<double> *>(mpp->mObjectives.at(0).get());
    int count = obj->mCounters.mFuncCalls;
    file << count;
}

bool testRosenbrockBench(std::shared_ptr<BM> bm, double eps, std::ofstream& file) {
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
    searchMethod.getOptions().mMaxStepsNumber = 100000;
    searchMethod.getOptions().mMinGrad = eps;
    searchMethod.getOptions().mHLB = eps;
    searchMethod.getOptions().mEps = eps;
    
    double x[dim];
    snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
    
    double v;
    
    searchMethod.search(x, v);
    
    file << v << "\t";
    file << std::abs(v - bm->getGlobMinY()) << "\t";
    
    bool isClose = std::abs(v - bm->getGlobMinY()) < eps;
    file << isClose << "\t";
    
    auto obj = dynamic_cast<COMPI::FuncCnt<double> *>(mpp->mObjectives.at(0).get());
    int count = obj->mCounters.mFuncCalls;
    file << count << "\t";
    
    
    // Retry without ortogonalization
    searchMethod.getOptions().mDoOrt = false;
    obj->mCounters.mFuncCalls = 0;
    snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
    searchMethod.search(x, v);
    isClose = std::abs(v - bm->getGlobMinY()) < eps;

    file << v << "\t";
    file << std::abs(v - bm->getGlobMinY()) << "\t";
    file << isClose << "\t";

    count = obj->mCounters.mFuncCalls;
    file << count;
}
