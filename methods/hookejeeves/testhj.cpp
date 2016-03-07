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

#include "hookjeeves.hpp"
#include "coorhjexplorer.hpp"
#include "rndhjexplorer.hpp"

class HJStopper : public LOCSEARCH::MHookeJeeves<double>::Stopper {
public:

    bool stopnow(double xdiff, double fdiff, double fval, int n) {
        //        std::cout << "n = " << n << "fdiff = " << fdiff << "\n";
        mCnt++;
        return false;
    }

    int mCnt = 0;
};

class QuadStopper : public LOCSEARCH::QuadLS<double>::Stopper {
public:

    bool stopnow(double s, int k, double vo, double vn) {
        mCnt++;
#if 0        
        std::cout << "s = " << s << ", k = " << k << "\n";
        std::cout << "vo = " << vo << ", vn = " << vn << "\n";
#endif
        if (s < 1e-3)
            return true;
        else if (k > 16)
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
    const int n = 10;
    OPTITEST::DejongProblemFactory fact(n, -4, 8);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);


    HJStopper stp;
#if 1    
    LOCSEARCH::CoorHJExplorer<double> explr(*mpp);
    explr.getOptions().mResetEveryTime = true;
#endif    
#if 0    
    LOCSEARCH::RndHJExplorer<double> explr(*mpp);
    explr.getOptions().mMaxTries = 12;
#endif

    QuadStopper lstp;
    LOCSEARCH::QuadLS<double> ls(*mpp, lstp);

#if 0
    LOCSEARCH::MHookeJeeves<double> desc(*mpp, stp, explr, nullptr);
#else    
    LOCSEARCH::MHookeJeeves<double> desc(*mpp, stp, explr, &ls);
#endif
    
    desc.getOptions().mInc = 1;
    desc.getOptions().mLambdaLB = desc.getOptions().mLambda;

    double x[n];

    for (int i = 0; i < n; i++)
        x[i] = i * 100 + 1;
    double v;
    bool rv = desc.search(x, v);
    std::cout << desc.about() << "\n";
    std::cout << "In " << stp.mCnt << " iterations found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
    SG_ASSERT(v <= 0.01);
    return 0;
}

