/* 
 * File:   tesgfsdesc.cpp
 * Author: medved
 *
 * Created on February 22, 2016, 3:17 PM
 */

#include <iostream>
#include <oneobj/contboxconstr/dejong.hpp>
#include <funccnt.hpp>
#include <box/boxutils.hpp>
#include <methods/lins/goldsec/goldsec.hpp>
#include "varcoordesc.hpp"

/*
 * 
 */
int main(int argc, char** argv) {
    const int n = 10;
    OPTITEST::DejongProblemFactory fact(n, -4, 8);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);


    LOCSEARCH::VarCoorDesc<double> desc(*mpp);
    desc.getOptions().mHInit = .1;
    desc.getOptions().mDoTracing = true;
    LOCSEARCH::GoldenSecLS<double>* locs = new LOCSEARCH::GoldenSecLS<double>(*mpp);
    locs->getOptions().mDoTracing = true;
    locs->getOptions().mSInit = 0.001;
    locs->getOptions().mMaxBackSteps = 1;
    desc.getLineSearch().reset(locs);

    double x[n];
    snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
    double v;
    int steps = 0;
    auto wtch = [&steps](double fval, const double* x, const std::vector<double>& gran, int stpn) {
        steps++;
    };
    desc.getWatchers().push_back(wtch);
    bool rv = desc.search(x, v);
    std::cout << desc.about() << "\n";
    std::cout << "In " << steps << " iterations found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
    SG_ASSERT(v <= 0.01);

    return 0;
}

