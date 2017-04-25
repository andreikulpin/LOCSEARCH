/* 
 * File:   tesgfsdesc.cpp
 * Author: medved
 *
 * Created on February 22, 2016, 3:17 PM
 */

#include <iostream>
#include <box/boxutils.hpp>
#include <oneobj/contboxconstr/dejong.hpp>
#include <oneobj/contboxconstr/rosenbrock.hpp>
#include <oneobj/contboxconstr/ackley1.hpp>
#include <funccnt.hpp>
#include <methods/lins/goldsec/goldsec.hpp>
#include <methods/lins/smartls/smartls.hpp>
#include "advancedcoordescent.hpp"

/*
 * 
 */
int main(int argc, char** argv) {
    const int n = 10;
    //OPTITEST::DejongProblemFactory fact(n, -4, 8);
    
    
    //OPTITEST::Ackley1ProblemFactory fact(std::vector<std::pair<double,double>>(n, std::pair<double, double>(-4,8)));
    OPTITEST::RosenbrockProblemFactory fact(n, -4, 8);
    COMPI::MPProblem<double> *mpp = fact.getProblem();
    COMPI::FuncCnt<double> *obj = new COMPI::FuncCnt<double>(*mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);
    LOCSEARCH::AdvancedCoordinateDescent<double> desc(*mpp);
#if 1    
    LOCSEARCH::GoldenSecLS<double>* locs = new LOCSEARCH::GoldenSecLS<double>(*mpp);
    locs->getOptions().mSInit = 0.1;
    locs->getOptions().mDelta = 0.02;
    locs->getOptions().mMaxBackSteps = 16;
    locs->getOptions().mDoTracing = true;
#endif
#if  0  
    LOCSEARCH::SmartLS<double> *locs = new LOCSEARCH::SmartLS<double>(*mpp);
    locs->getOptions().mSInit = 1e-4;
    locs->getOptions().mDoTracing = true;
    locs->getOptions().mMaxFailStepsBack = 0;
#endif    
    desc.getLineSearch().reset(locs);    
    desc.getOptions().mHInit = .1;
    desc.getOptions().mDoTracing = true;
    desc.getOptions().mSearchType = LOCSEARCH::AdvancedCoordinateDescent<double>::SearchTypes::PSEUDO_GRAD;
    //desc.getOptions().mSearchType = LOCSEARCH::AdvancedCoordinateDescent<double>::SearchTypes::HOOKE_JEEVES;
    //desc.getOptions().mSearchType = LOCSEARCH::AdvancedCoordinateDescent<double>::SearchTypes::NO_DESCENT;
    desc.getOptions().mVicinityAdaptation = LOCSEARCH::AdvancedCoordinateDescent<double>::VARIABLE_ADAPTATION;
    //desc.getOptions().mVicinityAdaptation = LOCSEARCH::AdvancedCoordinateDescent<double>::UNIFORM_ADAPTATION;
    double x[n];
    snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
    double v;
    bool rv = desc.search(x, v);
    std::cout << desc.about() << "\n";
    std::cout << "Found v = " << v << "\n";
    std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
    SG_ASSERT(v <= 0.01);

    return 0;
}

