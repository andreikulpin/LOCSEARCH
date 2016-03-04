/* 
 * File:   tesgfsdesc.cpp
 * Author: medved
 *
 * Created on February 22, 2016, 3:17 PM
 */

#include <iostream>
#include <oneobj/contboxconstr/dejong.hpp>
#include <funccnt.hpp>
#include <methods/lins/dichotls/dichotls.hpp>
#include <methods/lins/quadls/quadls.hpp>
#include "coordesc.hpp"

class CoorStopper : public LOCSEARCH::CoorDesc<double>::Stopper {
public:

    bool stopnow(double xdiff, double fdiff, double gran, double fval, int n) {
        mCnt++;
        return false;
    }

    int mCnt = 0;
};


/*
 * 
 */
int main(int argc, char** argv) {
    const int n = 1000;
    OPTITEST::DejongProblemFactory fact(n, -4, 8);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);


    CoorStopper stp;
    LOCSEARCH::CoorDesc<double> desc(*mpp, stp);

    desc.getOptions().mHInit = .1;

    double x[n];

    for (int i = 0; i < n; i++)
        x[i] = i;
    double v;
    bool rv = desc.search(x, v);
    std::cout <<  desc.about() << "\n";
    std::cout << "In " << stp.mCnt << " iterations found v = " << v << "\n";
    //std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";

    return 0;
}

