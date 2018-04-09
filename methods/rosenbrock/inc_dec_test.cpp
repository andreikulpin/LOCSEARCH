/* 
 * File:   inc_dec_test.cpp
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
    float minInc = 4.0;
    float maxInc = 7.0;
    float incStep = 0.1;
    
    float minDec = - 0.7;
    float maxDec = - 0.3;
    float decStep = 0.005;
    
    for (float inc = minInc; inc < maxInc; inc += incStep) {
        for (float dec = minDec; dec < maxDec; dec += decStep) {
            const int n = 50;

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
            desc.getOptions().mInc = inc;
            desc.getOptions().mDec = dec;
            desc.getOptions().maxUnsuccessStepsNumber = 50;
            desc.getOptions().maxStepsNumber = 100000;

            double x[n];
            snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);

            if (argc > 2) {
                desc.getOptions().mEps = atof(argv[2]);
            }

            double v;
            bool rv = desc.search(x, v);
            
            std::cout << "" << obj->mCounters.mFuncCalls << " ";
        }
        
        std::cout << "\n";
    }
    return 0;
}

