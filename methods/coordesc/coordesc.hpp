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
#include <solver.hpp>
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
    template <typename FT> class CoorDesc : public COMPI::Solver <FT> {
    public:

        /**
         * Determines stopping conditions
         * @param fval current best value found
         * @param x current best point
         * @param stpn step number
         * @return true if the search should stop
         */
        using Stopper = std::function<bool(FT fval, const FT* x, int stpn) >;

        /**
         * Watches the current step
         * @param fval current best value found
         * @param current best point
         * @param gran - current granularity 
         * @param stpn step number
         */
        using Watcher = std::function<void(FT fval, const FT* x, FT gran, int stpn) >;

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
         */
        CoorDesc(const COMPI::MPProblem<FT>& prob) :
        mProblem(prob) {
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
            FT h = mOptions.mHInit;
            const snowgoose::Box<double>& box = *(mProblem.mBox);

            // One step
            auto step = [&] () {
                int dir = 1;
                for (int i = 0; i < n;) {
                    FT y = x[i] + dir * h;
                    y = SGMAX(y, box.mA[i]);
                    y = SGMIN(y, box.mB[i]);
                    FT tmp = x[i];
                    x[i] = y;
                    FT fn = obj->func(x);
                    if (fn >= fcur) {
                        x[i] = tmp;
                        if (dir == 1)
                            dir = -1;
                        else {
                            dir = 1;
                            i++;
                        }
                    } else {
                        fcur = fn;
                        dir = 1;
                        i++;
                    }
                }
            };

            int sn = 0;
            bool br = false;
            while (!br) {
                sn++;
                FT fold = fcur;
                step();
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
                    w(fcur, x, h, sn);
                }

                for (auto s : mStoppers) {
                    if (s(fcur, x, sn)) {
                        br = true;
                        break;
                    }
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
            os << "Upper bound on the step size = " << mOptions.mHUB << "\n";
            os << "Lower bound on the step size = " << mOptions.mHLB << "\n";
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
         * Retrieve stoppers vector reference
         * @return stoppers vector reference
         */
        std::vector<Stopper>& getStoppers() {
            return mStoppers;
        }

        /**
         * Get watchers' vector
         * @return watchers vector
         */
        std::vector<Watcher>& getWatchers() {
            return mWatchers;
        }

    private:
        const COMPI::MPProblem<FT>& mProblem;
        Options mOptions;
        std::vector<Watcher> mWatchers;
        std::vector<Stopper> mStoppers;
    };
}

#endif /* GFSDESC_HPP */

