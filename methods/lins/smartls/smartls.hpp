/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   wolfels.hpp
 * Author: posypkin
 *
 * Created on April 21, 2017, 9:05 PM
 */

#ifndef SMARTLS_HPP
#define SMARTLS_HPP

#include <vector>
#include <limits>
#include <sstream>
#include <common/lineseach.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>


namespace LOCSEARCH {

    /**
     * Line search based on Wolfe's condition
     */
    template <class FT> class SmartLS : public LineSearch<FT> {
    public:

        struct Options {
            /**
             * Initial step
             */
            FT mSInit = 0.1;
            /**
             * Increase factor
             */
            FT mInc = 2;
            /**
             * Acceptance coefficient in Wolfe's condition (should be in (0,1))
             */
            FT mDec = 0.5;
            /**
             * Upper bound on unsuccessfull steps back
             */
            FT mMaxFailStepsBack = 1;
            /**
             * Trace on/off
             */
            bool mDoTracing = false;
            /**
             * Interactive search
             */
            bool mAskUser = false;

        };

        /**
         * The constructor
         *
         * @param prob optimization problem
         */
        SmartLS(const COMPI::MPProblem<FT>& prob) :
        mProblem(prob) {
            mXnew.resize(prob.mVarTypes.size());
            mDir.resize(prob.mVarTypes.size());
        }

        bool search(const FT* d, FT* x, FT& v) override {
            FT l = mOptions.mSInit;
            FT lbest = 0;
            FT fbest = v;
            FT lbase = mOptions.mSInit;
            const int n = mProblem.mVarTypes.size();
            FT* xnew = mXnew.data();
            const snowgoose::Box<FT>& box = *(mProblem.mBox);
            auto obj = mProblem.mObjectives.at(0);
            bool br = false;
            int sinc = 0;
            int sdec = 0;
            int fdec = 0;
            while (!br) {
                for (int j = 0; j < n; j++) {
                    FT nalpha = l;
                    if (d[j] > 0) {
                        nalpha = (box.mB[j] - x[j]) / d[j];
                    } else if (d[j] < 0) {
                        nalpha = (box.mA[j] - x[j]) / d[j];
                    }
                    l = SGMIN(l, nalpha);
                }
                snowgoose::VecUtils::vecSaxpy(n, x, d, l, xnew);

                FT fn = obj->func(xnew);
                if (fn < fbest) {
                    if (mOptions.mDoTracing)
                        std::cout << "SUCCESS: " << " l = " << l << ", fn = " << fn << ", fbest (old)  = " << fbest << "\n";
                    fbest = fn;
                    lbest = l;
                    if (fdec == 0) {
                        l *= mOptions.mInc;
                        sinc ++;
                    } else {
                        l *= mOptions.mDec;
                        sdec ++;
                    }
                } else {
                    if (mOptions.mDoTracing)
                        std::cout << "FAIL: " << " l = " << l << ", fn = " << fn << ", fbest = " << fbest << "\n";
                    if (sinc == 0) {
                        if((sdec == 0) && (fdec < mOptions.mMaxFailStepsBack)) {
                            l *= mOptions.mDec;
                            fdec ++;
                        } else
                            br = true;
                    } else
                        br = true;
                }

                if (br) {
                    if (mOptions.mAskUser) {
                        FT nl = askForL();
                        if (nl == 0)
                            break;
                        else {
                            l = nl;
                            br = false;
                        }
                    } else
                        break;
                }
            }


            if (fbest < v) {
                snowgoose::VecUtils::vecSaxpy(n, x, d, lbest, x);
                v = fbest;
                return true;
            } else
                return false;
        }

        std::string about() const {
            std::ostringstream os;
            os << "Smart Line Search with \n";
            
            os << "Initial step = " << mOptions.mSInit << "\n";
            os << "Increment multiplier = " << mOptions.mInc << "\n";
            os << "Decrement multiplier = " << mOptions.mDec << "\n";
            os << "Maximal back steps = " << mOptions.mMaxFailStepsBack << "\n";
            return os.str();
        }

        /**
         * Retrieve options reference
         * @return options reference
         */
        Options& getOptions() {
            return mOptions;
        }

    private:

        FT askForL() const {
            FT l;
            std::cout << "Suggest l: ";
            std::cin >> l;
            return l;
        }
        std::vector<FT> mXnew;
        std::vector<FT> mDir;
        const COMPI::MPProblem<FT>& mProblem;
        Options mOptions;
    };
}


#endif /* SMARTLS_HPP */

