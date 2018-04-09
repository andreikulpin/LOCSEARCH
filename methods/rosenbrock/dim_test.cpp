/* 
 * File:   dim_test.cpp
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
#include "rosenbrock_method.hpp"

/*
 * 
 */
int main(int argc, char** argv) {
    int minN = 2;
    int maxN = 100;
    int stepN = 1;
    
    for (int i = minN; i <= maxN; i += stepN) {
        const int n = i;

        OPTITEST::RosenbrockProblemFactory fact(n, -2, 5);
        COMPI::MPProblem<double> *mpp = fact.getProblem();
        auto obj = std::make_shared<COMPI::FuncCnt<double>>(mpp->mObjectives.at(0));
        mpp->mObjectives.pop_back();
        mpp->mObjectives.push_back(obj);

        double initH[n];
        snowgoose::VecUtils::vecSet(n, .1, initH);

        LOCSEARCH::RosenbrockMethod<double> desc(*mpp);
        desc.getOptions().mHInit = initH;
        desc.getOptions().mDoTracing = false;
        desc.getOptions().mEps = 1e-10;
        desc.getOptions().mInc = 5.0;
        desc.getOptions().mDec = - 0.5;
        desc.getOptions().maxUnsuccessStepsNumber = 50;
        desc.getOptions().maxStepsNumber = 100000;

        double x[n];
        snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);

        double v;
        bool rv = desc.search(x, v);

        std::cout << "" << n << " ";
        std::cout << "" << v << " ";
        std::cout << "" << obj->mCounters.mFuncCalls << "\n";
    }

    return 0;
}

