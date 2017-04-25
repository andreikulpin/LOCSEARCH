/* 
 * File:   tesgfsdesc.cpp
 * Author: medved
 *
 * Created on February 22, 2016, 3:17 PM
 */

#include <iostream>
#include <oneobj/contboxconstr/dejong.hpp>
#include <funccnt.hpp>
#include "smartls.hpp"



/*
 * 
 */
int main(int argc, char** argv) {
    const int n = 2;
    OPTITEST::DejongProblemFactory fact(n, -4, 4);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();    
    mpp->mObjectives.push_back(obj);
    
    LOCSEARCH::SmartLS<double> ls(*mpp);
    ls.getOptions().mDoTracing = true;

    double x[n], d[n];

    for (int i = 0; i < n; i++) {
        x[i] = 1;
        d[i] = -1;
    }
    double v = obj->func(x);
    bool rv = ls.search(d, x, v);
    std::cout << "Found v = " << v << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";

    return 0;
}

