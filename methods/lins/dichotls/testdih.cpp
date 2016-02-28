/* 
 * File:   tesgfsdesc.cpp
 * Author: medved
 *
 * Created on February 22, 2016, 3:17 PM
 */

#include <iostream>
#include <oneobj/contboxconstr/dejong.hpp>
#include <funccnt.hpp>
#include "dichotls.hpp"


class MyStopper : public LOCSEARCH::DichotLS<double>::Stopper {
public:

    bool stopnow(double s, int k, double vo, double vn){
        mCnt++;
        std::cout << "k = " << k << ", s = " << s << ", vo = " << vo << ", vn = " << vn << "\n";
        if (s < 1e-5)
            return true;
        else if(k > 128)
            return true;
        else
            return false;
    }

    int mCnt = 0;
};

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
    
    MyStopper stp;
    LOCSEARCH::DichotLS<double> ls(*mpp, stp);

    //desc.getOptions().mOnlyCoordinateDescent = true;

    double x[n], d[n];

    for (int i = 0; i < n; i++) {
        x[i] = 1;
        d[i] = -1;
    }
    double v;
    bool rv = ls.search(d, x, v);
    std::cout << "In " << stp.mCnt << " iterations found v = " << v << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";

    return 0;
}

