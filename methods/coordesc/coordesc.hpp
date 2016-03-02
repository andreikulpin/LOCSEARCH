/* 
 * File:   bbboxdesc.hpp
 * Author: medved
 *
 * Created on November 3, 2015, 5:05 PM
 */

#ifndef COORDESC_HPP
#define  COORDESC_HPP

#include <common/locsearch.hpp>
#include <common/lineseach.hpp>
#include <common/dummyls.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include <common/sgerrcheck.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>

namespace LOCSEARCH {

    /**
     * Simple coordinate descen for box constrained problems
     */
    template <typename FT> class CoorDesc : public LocalSearch <FT> {
    public:

        /**
         * Determines stopping conditions
         */
        class Stopper {
        public:
            /**
             * Returns true when the search should stop
             * @param xdiff difference between old and new x
             * @param fdiff difference between old and new f value
             * @param gran current granularity
             * @param fval function value
             * @param n current step number 
             */
            virtual bool stopnow(FT xdiff, FT fdiff, FT gran, FT fval, int n) = 0;
        };

        /**
         * Options for Gradient Box Descent method
         */
        struct Options {
            /**
             * Initial value of granularity
             */
            FT mHInit = 0.01;

            /**
             * Increase in the case of success
             */
            FT mInc = 1.75;
            /**
             * Decrease in the case of failure
             */
            FT mDec = 0.5;
            /**
             * Lower bound for granularity
             */
            FT mHLB = 1e-08;
            /**
             * Upper bound on granularity
             */
            FT mHUB = 1e+02;
        };

        /**
         * The constructor
         * @param prob - reference to the problem
         * @param stopper - reference to the stopper
         * @param ls - pointer to the line search
         */
        CoorDesc(const COMPI::MPProblem<FT>& prob, Stopper& stopper) :
        mProblem(prob),
        mStopper(stopper) {
            unsigned int typ = COMPI::MPUtils::getProblemType(prob);
            SG_ASSERT(typ = COMPI::MPUtils::ProblemTypes::BOXCONSTR | COMPI::MPUtils::ProblemTypes::CONTINUOUS | COMPI::MPUtils::ProblemTypes::SINGLEOBJ);
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
            FT h = mOptions.mHInit;
            const snowgoose::Box<double>& box = *(mProblem.mBox);

            // One step
            auto step = [&] () {
                bool succ = false;
                FT xd = 0;
                for (int i = 0; i < n; i++) {
                    FT y = x[i] - h;
                    if (y < box.mA[i]) {
                        y = box.mA[i];
                    }
                    FT tmp = x[i];
                    x[i] = y;
                    FT fn = obj->func(x);
                    if (fn >= fcur) {
                        x[i] = tmp;
                    } else {
                        fcur = fn;
                        xd += h * h;
                        continue;
                    }

                    y = x[i] + h;
                    if (y > box.mB[i]) {
                        y = box.mB[i];
                    }
                    tmp = x[i];
                    x[i] = y;
                    fn = obj->func(x);
                    if (fn >= fcur) {
                        x[i] = tmp;
                    } else {
                        fcur = fn;
                        xd += h * h;
                        continue;
                    }
                }
                return xd;
            };

            int sn = 0;
            for (;;) {
                sn++;
                FT fold = fcur;
                FT xdiff = step();
                FT fdiff = fold - fcur;
                if (fcur < fold) {
                    rv = true;
                    h *= mOptions.mInc;
                    h = SGMIN(h, mOptions.mHUB);
                } else {
                    if (h <= mOptions.mHLB)
                        break;
                    h *= mOptions.mDec;
                }
                if (mStopper.stopnow(xdiff, fdiff, h, fcur, sn)) {
                    break;
                }

            }
            v = fcur;
            return rv;
        }

        /**
         * Retrieve options
         * @return options
         */
        Options & getOptions() {
            return mOptions;
        }

    private:

        /*
        void project(FT* x) {
            int n = LocalOptimizer<FT>::getObjective()->getDim();
            for (int i = 0; i < n; i++) {
                if (x[i] < LocalBoxOptimizer<FT>::mBox.mA[i])
                    x[i] = LocalBoxOptimizer<FT>::mBox.mA[i];
                else if (x[i] > LocalBoxOptimizer<FT>::mBox.mB[i])
                    x[i] = LocalBoxOptimizer<FT>::mBox.mB[i];
            }
        }
         */

        Stopper &mStopper;
        const COMPI::MPProblem<FT>& mProblem;
        LineSearch<FT>* mLS;
        Options mOptions;
    };
}

#endif /* GFSDESC_HPP */

