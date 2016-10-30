/*
 * Coordinate descent for constrained problems that uses projection
 */

/* 
 * File:   conscoordesc.hpp
 * Author: posypkin
 *
 * Created on October 28, 2016, 8:43 PM
 */

#ifndef CONSCOORDESC_HPP
#define CONSCOORDESC_HPP
#include <sstream>
#include <functional>
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
     * Simple coordinate descen for constrained problems
     */
    template <typename FT> class ProjCoorDesc : public LocalSearch <FT> {
    public:

        /**
         * Determines stopping conditions
         */
        typedef std::function<bool(FT xdiff, FT fdiff, FT gran, FT fval, int n)> Stopper;
        

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
         * @param projector - projector function
         */
        ProjCoorDesc(const COMPI::MPProblem<FT>& prob,
                const Stopper& stopper,
                const std::function<void(FT*)>& projector) :
        mProblem(prob),
        mStopper(stopper),
        mProjector(projector) {
            unsigned int typ = COMPI::MPUtils::getProblemType(prob);
            SG_ASSERT(typ == COMPI::MPUtils::ProblemTypes::CONTINUOUS | COMPI::MPUtils::ProblemTypes::SINGLEOBJ);
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
            mProjector(x);
            FT fcur = obj->func(x);
            int n = mProblem.mVarTypes.size();
            FT h = mOptions.mHInit;
            const snowgoose::Box<double>& box = *(mProblem.mBox);
            FT y[n];
            FT z[n];

            // One step
            auto step = [&] () {
                snowgoose::VecUtils::vecCopy(n, x, z);
                for (int i = 0; i < n; i++) {
                    snowgoose::VecUtils::vecCopy(n, x, y);
                    x[i] -= h;
                    mProjector(x);
                    FT fn = obj->func(x);
                    if (fn >= fcur) {
                        snowgoose::VecUtils::vecCopy(n, y, x);
                    } else {
                        fcur = fn;
                        continue;
                    }

                    x[i] += h;
                    mProjector(x);
                    fn = obj->func(x);
                    if (fn >= fcur) {
                        snowgoose::VecUtils::vecCopy(n, y, x);
                    } else {
                        fcur = fn;                        
                    }
                }
                FT xd = snowgoose::VecUtils::vecDist(n, z, x);
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
                if (mStopper(xdiff, fdiff, h, fcur, sn)) {
                    break;
                }

            }
            v = fcur;
            return rv;
        }

        std::string about() const {
            std::ostringstream os;
            os << "Projected coordinate descent method\n";
            os << "Initial step = " << mOptions.mHInit << "\n";
            os << "Increment multiplier = " << mOptions.mInc << "\n";
            os << "Decrement multiplier = " << mOptions.mDec << "\n";
            os << "Upper bound on the step = " << mOptions.mHUB << "\n";
            os << "Lower bound on the step = " << mOptions.mHLB << "\n";
            return os.str();
        }

        /**
         * Retrieve options
         * @return options
         */
        Options & getOptions() {
            return mOptions;
        }

    private:

        Stopper mStopper;
        const COMPI::MPProblem<FT>& mProblem;
        Options mOptions;
        std::function<void(FT* x)> mProjector;
    };
}




#endif /* CONSCOORDESC_HPP */

