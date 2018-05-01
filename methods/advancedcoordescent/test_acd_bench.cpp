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
#include "advancedcoordescent.hpp"
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
    //locs->getOptions().mDoTracing = true;
    
    LOCSEARCH::AdvancedCoordinateDescent<double> searchMethod(*mpp);
    searchMethod.getLineSearch().reset(locs);    
    searchMethod.getOptions().mHInit = .1;
    //searchMethod.getOptions().mDoTracing = true;
    searchMethod.getOptions().mGradLB = 0;
    
    double x[dim];
    for (int i = 0; i < dim; i++) {
        double a = bm->getBounds()[i].first;
        double b = bm->getBounds()[i].second;
        x[i] = (b + a) / 2.0;
    }
    
    double v;
    
    std::cout << "*************Testing benchmark**********" << std::endl;
    searchMethod.search(x, v);
    std::cout << "v = " << v << std::endl;
    std::cout << "the difference is " << v - bm->getGlobMinY() << std::endl;
}

int main(int argc, char** argv) {
    const int dim = argc > 1 ? atoi(argv[1]) : 50;
    double eps = argc > 2 ? atof(argv[2]) : 0.01;
    
    Benchmarks<double> tests;
    for (auto bm : tests) {
        testBench(bm, eps);
    }
    return 0;
}

