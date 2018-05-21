/* 
 * File:   testrosenbrock3.cpp
 * Author: kulpin
 */

#include <iostream>
#include <box/boxutils.hpp>
#include <oneobj/contboxconstr/dejong.hpp>
#include <oneobj/contboxconstr/rosenbrock.hpp>
#include <oneobj/contboxconstr/ackley1.hpp>
#include <funccnt.hpp>
#include <methods/lins/goldsec/goldsec.hpp>
#include <methods/lins/smartls/smartls.hpp>
#include "rosenbrockmethod.hpp"

/*
 * 
 */
int main(int argc, char** argv) {
    const int dim = argc > 1 ? atoi(argv[1]) : 50;
    const double eps = argc > 2 ? atof(argv[2]) : 0.01;
    
    OPTITEST::RosenbrockProblemFactory fact(dim, -30, 30);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    auto obj = std::make_shared<COMPI::FuncCnt<double>>(mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);
    
    double initH[dim];
    snowgoose::VecUtils::vecSet(dim, .1, initH);
    
    LOCSEARCH::RosenbrockMethod<double> searchMethod(*mpp);
    searchMethod.getOptions().mHInit = std::vector<double>(initH, initH + dim);
    searchMethod.getOptions().mDoTracing = false;
    searchMethod.getOptions().mInc = 5.0;
    searchMethod.getOptions().mDec = 0.5;
    searchMethod.getOptions().mMaxStepsNumber = 100000;
    searchMethod.getOptions().mMinGrad = eps;
    searchMethod.getOptions().mHLB = searchMethod.getOptions().mMinGrad;
    
    double x[dim];
    snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
    
    double v;
    bool rv = searchMethod.search(x, v);
    
    std::cout << searchMethod.about() << "\n";
    std::cout << "Found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(dim, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
    SG_ASSERT(v <= 0.01);

    return 0;
}

