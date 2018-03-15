/* 
 * File:   testrosenbrock.cpp
 * Author: Kulpin
 */

#include <iostream>
#include <box/boxutils.hpp>
#include <oneobj/contboxconstr/dejong.hpp>
#include <oneobj/contboxconstr/rosenbrock.hpp>
#include <oneobj/contboxconstr/ackley1.hpp>
#include <funccnt.hpp>
#include <methods/lins/goldsec/goldsec.hpp>
#include <methods/lins/smartls/smartls.hpp>
#include "rosenbrock_method.hpp"
#include "testfunc.hpp"

/*
 * 
 */
int main(int argc, char** argv) {
    const int n = 2;
    
    OPTITEST::TestProblemFactory fact(n, -4, 8);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    auto obj = std::make_shared<COMPI::FuncCnt<double>>(mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);

    double initH[] = {1., 2.};
    
    LOCSEARCH::RosenbrockMethod<double> desc(*mpp);   
    desc.getOptions().mHInit = initH;
    desc.getOptions().mDoTracing = true;
    desc.getOptions().mEps = 0.6;
    desc.getOptions().mInc = 2;
    desc.getOptions().mDec = - 0.5;

    double x[n];
    x[0] = 8.;
    x[1] = 9.;

    double v;
    bool rv = desc.search(x, v);

    std::cout << desc.about() << "\n";
    std::cout << "Found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
    return 0;
}

