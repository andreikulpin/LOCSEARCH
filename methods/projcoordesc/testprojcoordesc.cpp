/* 
 * File:   tesgfsdesc.cpp
 * Author: medved
 *
 * Created on February 22, 2016, 3:17 PM
 */

#include <iostream>
#include <cmath>
#include <oneobj/contboxconstr/dejong.hpp>
#include <funccnt.hpp>
#include <methods/lins/dichotls/dichotls.hpp>
#include <methods/lins/quadls/quadls.hpp>
#include "projcoordesc.hpp"



/*
 * 
 */
int main(int argc, char** argv) {
    const int n = 10;
    OPTITEST::DejongProblemFactory fact(n, 4, 8);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);

    const snowgoose::Box<double>& box = *(mpp->mBox);
    int cnt = 0;
    auto stopper = [&](double xdiff, double fdiff, double gran, double fval, int n){cnt ++; return false;};
    auto projector = [&box] (double* x) {
        snowgoose::BoxUtils::project(x, box);
    };
    
    LOCSEARCH::ProjCoorDesc<double> desc(*mpp, stopper, projector);

    desc.getOptions().mHInit = .1;

    double x[n];

    for (int i = 0; i < n; i++)
        x[i] = i * 100 + 1;
    double v;
    bool rv = desc.search(x, v);
    std::cout << desc.about() << "\n";
    std::cout << "In " << cnt << " iterations found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
    SG_ASSERT(std::fabs(v - 160) <= 0.01);
    
    return 0;
}

