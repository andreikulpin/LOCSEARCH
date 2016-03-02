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
#include "gfsdesc.hpp"

class GFSStopper : public LOCSEARCH::GFSDesc<double>::Stopper {
public:

    bool stopnow(double xdiff, double fdiff, double gmin, double fval, int n) {
        mCnt++;
        if (gmin < 1e-3)
            return true;
        else
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

class DichStopper : public LOCSEARCH::DichotLS<double>::Stopper {
public:

    bool stopnow(double s, int k, double vo, double vn) {
        mCnt++;
#if 0        
        std::cout << "s = " << s << ", k = " << k << "\n";
        std::cout << "vo = " << vo << ", vn = " << vn << "\n";
#endif
        if (s < 1e-3)
            return true;
        else if (k > 512)
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
    const int n = 1000;
    OPTITEST::DejongProblemFactory fact(n, -4, 8);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);

#if 0    
    DichStopper lstp;
    LOCSEARCH::DichotLS<double> ls(*mpp, lstp);
    ls.getOptions().mSInit = .1;
    ls.getOptions().mAccelerate = 1.1;
    ls.getOptions().mSlowDown = 0.5;
#else
    QuadStopper lstp;
    LOCSEARCH::QuadLS<double> ls(*mpp, lstp);

#endif    

    GFSStopper stp;

#if 1      
    LOCSEARCH::GFSDesc<double> desc(*mpp, stp, &ls);
#else
    LOCSEARCH::GFSDesc<double> desc(*mpp, stp);
#endif    

#if 0
    desc.getOptions().mOnlyCoordinateDescent = true;
#endif

    desc.getOptions().mHInit = .1;

    double x[n];

    for (int i = 0; i < n; i++)
        x[i] = i;
    double v;
    bool rv = desc.search(x, v);
    std::cout << "In " << stp.mCnt << " iterations found v = " << v << "\n";
    //std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";

    return 0;
}

