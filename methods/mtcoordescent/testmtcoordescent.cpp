/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   testmtcoordescent.cpp
 * Author: mposypkin
 *
 * Created on October 26, 2017, 12:30 PM
 */

#include <iostream>
#include <chrono>
#include <box/boxutils.hpp>
#include <oneobj/contboxconstr/dejong.hpp>
#include <oneobj/contboxconstr/rosenbrock.hpp>
#include <oneobj/contboxconstr/ackley1.hpp>
#include <funccnt.hpp>
#include <methods/lins/goldsec/goldsec.hpp>
#include <methods/lins/smartls/smartls.hpp>
#include "ctcoordescent.hpp"

/*
 * 
 */
int main(int argc, char** argv) {
    const int n = 500;
    //OPTITEST::DejongProblemFactory fact(n, -4, 8);


    //OPTITEST::Ackley1ProblemFactory fact(std::vector<std::pair<double,double>>(n, std::pair<double, double>(-4,8)));
    OPTITEST::RosenbrockProblemFactory fact(n, -2, 5);
    COMPI::MPProblem<double> *mpp = fact.getProblem();

    /*
    auto obj = std::make_shared<COMPI::FuncCnt<double>>(mpp->mObjectives.at(0));
    mpp->mObjectives.pop_back();
    mpp->mObjectives.push_back(obj);
     */

    LOCSEARCH::CTCoordinateDescent<double> desc(*mpp);
    desc.getOptions().mHInit = .1;
    //desc.getOptions().mHLB = 1e-10;
    //desc.getOptions().mParallelMode = false;
    desc.getOptions().mNumThreads = 6;
    desc.getOptions().mGradLB = 0;

    double x[n];
    snowgoose::BoxUtils::getCenter(*(mpp->mBox), x);
    double v;
    auto start = std::chrono::steady_clock::now();
    bool rv = desc.search(x, v);
    auto end = std::chrono::steady_clock::now();
    auto diff = end - start;
    std::cout << desc.about() << "\n";
    std::cout << "Found v = " << mpp->mObjectives.at(0)->func(x) << "\n";
    //std::cout << " at " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
    //std::cout << "Number of objective calls is " << obj->mCounters.mFuncCalls << "\n";
    //SG_ASSERT(v <= 0.01);
    std::cout << "Execution time = " << std::chrono::duration <double, std::milli> (diff).count() << " ms" << std::endl;
    return 0;
}



