/* 
 * File:   bbboxdesc.hpp
 * Author: medved
 *
 * Created on November 3, 2015, 5:05 PM
 */

#ifndef COORDESC_HPP
#define  COORDESC_HPP

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
     * Simple coordinate descent for box constrained problems
     */
    template <typename FT> class CoorDesc : public LocalSearch <FT> {
    public:

        /**
         * Determines stopping conditions
         * @param xdiff - distance between next and previous x
         * @param fdiff - difference between next and previous f value
         * @param gran - current granularity
         * @param n - current step number
         * @return true if search should stop
         */
        typedef std::function<bool(FT xdiff, FT fdiff, FT gran, FT fval, int n) > Stopper;

        /**
         * Watches the step data
         * @param xdiff - distance between next and previous x
         * @param fdiff - difference between next and previous f value
         * @param gran - current granularity
         * @param n - current step number
         */
        typedef std::function<void(FT xdiff, FT fdiff, FT gran, FT fval, int n) > Watcher;

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
            /**
             * Shifts along coordinate directions, by default it is a vector of ones. Used if we need different shifts 
             * along different coordinates.
             */
            std::vector<FT> mShifts;
        };

        /**
         * The constructor
         * @param prob - reference to the problem
         * @param stopper - reference to the stopper (copies stopper)
         * @param ls - pointer to the line search
         */
        CoorDesc(const COMPI::MPProblem<FT>& prob, const Stopper& stopper) :
        mProblem(prob),
        mStopper(stopper) {
            unsigned int typ = COMPI::MPUtils::getProblemType(prob);
            mOptions.mShifts.assign(prob.mVarTypes.size(), 1);
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
            FT h = mOptions.mHInit;
            const snowgoose::Box<double>& box = *(mProblem.mBox);

            // One step
            auto step = [&] () {
                FT xd = 0;
                for (int i = 0; i < n; i++) {
                    FT dh = h * mOptions.mShifts[i];
                    FT y = x[i] - dh;
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
                        xd += dh * dh;
                        continue;
                    }

                    y = x[i] + dh;
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
                        xd += dh * dh;
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

                for (auto w : mWatchers) {
                    w(xdiff, fdiff, h, fcur, sn);
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
            os << "Coordinate descent method\n";
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

        /**
         * Get watchers' vector
         * @return watchers vector
         */
        std::vector<Watcher>& getWatchers() {
            return mWatchers;
        }

    private:


        Stopper mStopper;
        const COMPI::MPProblem<FT>& mProblem;
        Options mOptions;
        std::vector<Watcher> mWatchers;

    };
}

#endif /* GFSDESC_HPP */

