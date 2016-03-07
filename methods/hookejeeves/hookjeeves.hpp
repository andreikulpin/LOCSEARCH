/* 
 * File:   bbboxdesc.hpp
 * Author: medved
 *
 * Created on November 3, 2015, 5:05 PM
 */

#ifndef HOOKEJEEVES_HPP
#define  HOOKEJEEVES_HPP

#include <sstream>
#include <common/locsearch.hpp>
#include <common/lineseach.hpp>
#include <common/dummyls.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include <common/sgerrcheck.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>

#include "hjexplorer.hpp"


namespace LOCSEARCH {

    /**
     * Modified Hooke-Jeeves Method
     */
    template <typename FT> class MHookeJeeves : public LocalSearch <FT> {
    public:

        struct Options {
            /**
             * Speedup multiplier
             */
            FT mLambda = 1;
            /**
             * Increment multiplier
             */
            FT mInc = 1;
            /**
             * Decrement multiplier
             */
            FT mDec = 0.5;
            /**
             * Lower bound on lambda
             */
            FT mLambdaLB = 0.01;
        };

        /**
         * Determines stopping conditions
         */
        class Stopper {
        public:
            /**
             * Returns true when the search should stop
             * @param xdiff difference between old and new x
             * @param fdiff difference between old and new f value
             * @param fval function value best at the moment
             * @param n current (big) step number 
             */
            virtual bool stopnow(FT xdiff, FT fdiff, FT fval, int n) = 0;
        };

        /**
         * The constructor
         * @param prob - reference to the problem
         * @param stopper - reference to the stopper
         * @param explorer - reference to the explorer
         * @param ls - pointer to the line search
         */
        MHookeJeeves(const COMPI::MPProblem<FT>& prob, Stopper& stopper, HJExplorer<FT>& explorer, LineSearch<FT>* ls = nullptr) :
        mProblem(prob),
        mStopper(stopper),
        mExplorer(explorer),
        mLS(ls) {
            unsigned int typ = COMPI::MPUtils::getProblemType(prob);
            SG_ASSERT(typ == COMPI::MPUtils::ProblemTypes::BOXCONSTR | COMPI::MPUtils::ProblemTypes::CONTINUOUS | COMPI::MPUtils::ProblemTypes::SINGLEOBJ);
        }

        /**
         * Perform search
         * @param x start point and result
         * @param v  the resulting value
         * @return true if search converged and false otherwise
         */
        bool search(FT* x, FT& v) {
            bool rv = false;

            COMPI::Functor<FT>* obj = mProblem.mObjectives.at(0);
            snowgoose::BoxUtils::project(x, *(mProblem.mBox));
            FT fcur = obj->func(x);
            int n = mProblem.mVarTypes.size();
            const snowgoose::Box<double>& box = *(mProblem.mBox);
            int sn = 0;
            FT lam = mOptions.mLambda;
            FT dir[n];
            FT xold[n];
            FT y[n];
            
            auto step = [&] (FT* x1, FT* x2, FT * x3) {
                if (mLS == nullptr) {
                    for (int i = 0; i < n; i++) {
                        x3[i] = x2[i] + lam * (x2[i] - x1[i]);
                    }
                } else {

                    FT vv;
                    snowgoose::VecUtils::vecSaxpy(n, x2, x1, -1., dir);
                    snowgoose::VecUtils::vecCopy(n, x2, x3);
                    mLS->search(dir, x3, vv);
                }
            };


            snowgoose::VecUtils::vecCopy(n, x, y);

            for (;;) {
                sn++;
                FT fnew = mExplorer.explore(y);
                if (fnew < fcur) {
                    rv = true;

                    FT fdiff = fcur - fnew;
                    fcur = fnew;
                    snowgoose::VecUtils::vecCopy(n, x, xold);
                    snowgoose::VecUtils::vecCopy(n, y, x);
                    FT xdiff = snowgoose::VecUtils::vecDist(n, xold, x);
                    step(xold, x, y);
                    lam *= mOptions.mInc;
                    if (mStopper.stopnow(xdiff, fdiff, fcur, sn))
                        break;

                } else {
                    if (lam > mOptions.mLambdaLB) {
                        lam *= mOptions.mDec;
                        snowgoose::VecUtils::vecCopy(n, x, y);
                    } else
                        break;
                }
            }
            v = fcur;
            return rv;
        }

        std::string about() const {
            std::ostringstream os;
            os << "Modified Hooke-Jeeves method\n";
            os << "Initial speedup coefficient: " << mOptions.mLambda << "\n";
            os << "Lower bound on speedup coefficient: " << mOptions.mLambdaLB << "\n";
            os << "Decrement multiplier: " << mOptions.mDec << "\n";
            os << "Increment multiplier: " << mOptions.mInc << "\n";
            os << "Explorer: \n" << mExplorer.about() << "\n";
            if (mLS != nullptr)
                os << "Line search: " << mLS->about() << "\n";
            return os.str();
        }

        Options& getOptions() {
            return mOptions;
        }

    private:


        Stopper &mStopper;
        HJExplorer<FT> &mExplorer;
        const COMPI::MPProblem<FT>& mProblem;
        Options mOptions;
        LineSearch<FT>* mLS;

    };
}

#endif 

