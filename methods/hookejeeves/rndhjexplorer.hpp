/* 
 * File:   hjexplorer.hpp
 * Author: medved
 *
 * Created on March 5, 2016, 3:23 PM
 */

#ifndef RNDHJEXPLORER_HPP
#define	RNDHJEXPLORER_HPP

#include <stdlib.h>
#include <mpproblem.hpp>
#include "hjexplorer.hpp"

namespace LOCSEARCH {

    /**
     * Random  exploration to be used in Hooke and Jeeves Method
     */
    template <class FT> class RndHJExplorer : public HJExplorer <FT> {
    public:

        struct Options {
            /**
             * Initial value of granularity
             */
            FT mHInit = 1;

            /**
             * Increase in the case of success
             */
            FT mInc = 1;
            /**
             * Decrease in the case of failure
             */
            FT mDec = 0.5;
            /**
             * Lower bound for granularity
             */
            FT mHLB = 1e-01;
            /**
             * Upper bound on granularity
             */
            FT mHUB = 1e+02;
            /**
             * maximal number of tries before improvement
             */
            unsigned int mMaxTries = 16;
            /**
             * Reset parameter in each invocation
             */
            bool mResetEveryTime = false;

        };

        RndHJExplorer(const COMPI::MPProblem<FT>& prob) : mProb(prob) {
            mH = mOptions.mHInit;
        }

        /**
         * Explore the vicinity of the 'x' point to find better value
         * @param x start vector on entry, resulting vector on exit 
         * @return new obtained value 
         */
        FT explore(FT* x) {
            COMPI::Functor<FT>* obj = mProb.mObjectives.at(0);
            int n = mProb.mVarTypes.size();
            const snowgoose::Box<double>& box = *(mProb.mBox);
            FT fcur = obj->func(x);
            FT xn[n];
            snowgoose::VecUtils::vecCopy(n, x, xn);
            reset();
            unsigned int cnt = 0;
            for (;;) {
                for (int i = 0; i < n; i++) {
                    FT r = 0.5 - (FT) rand() / (FT) RAND_MAX;
                    FT y = x[i] + mH * r;
                    if (y < box.mA[i]) {
                        y = box.mA[i];
                    }
                    if (y > box.mB[i]) {
                        y = box.mB[i];
                    }
                    xn[i] = y;
                }
                FT fn = obj->func(xn);
                if (fn < fcur) {
                    snowgoose::VecUtils::vecCopy(n, xn, x);
                    fcur = fn;
                    cnt = 0;
                    mH *= mOptions.mInc;
                    mH = SGMIN(mOptions.mHUB, mH);
                } else {
                    cnt++;
                    if (cnt > mOptions.mMaxTries) {
                        cnt = 0;
                        mH *= mOptions.mDec;
                        if (mH <= mOptions.mHLB)
                            break;
                    }
                }
                //std::cout << "cnt = " << cnt << ", h = " << mH << ", fcur = " << fcur << "\n";
            }
            return fcur;
        }

        std::string about() const {
            std::ostringstream os;
            os << "Randomized exporer\n";
            os << "Initial granularity: " << mOptions.mHInit << "\n";
            os << "Increment coefficient: " << mOptions.mInc << "\n";
            os << "Decrement coefficient: " << mOptions.mDec << "\n";
            os << "Lower bound on granularity: " << mOptions.mHLB << "\n";
            os << "Upper bound on granularity: " << mOptions.mHUB << "\n";
            os << "Max Tries: " << mOptions.mMaxTries << "\n";
            if (mOptions.mResetEveryTime)
                os << "Reset on every iteration\n";
            return os.str();
        }

        /**
         * Retrieve options
         * @return reference to options
         */
        Options& getOptions() {
            return mOptions;
        }

        /**
         * Reset the current vicinity size
         */
        void reset() {
            mH = mOptions.mHInit;
        }

    private:

        const COMPI::MPProblem<FT>& mProb;
        Options mOptions;
        FT mH;

    };

};


#endif	/* HJEXPLORER_HPP */

