/* 
 * File:   bbboxdesc.hpp
 * Author: medved
 *
 * Created on November 3, 2015, 5:05 PM
 */

#ifndef VARCOORDESC_HPP
#define  VARCOORDESC_HPP

#include <sstream>
#include <vector>
#include <functional>
#include <solver.hpp>
#include <common/lineseach.hpp>
#include <common/dummyls.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include <common/sgerrcheck.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>


namespace LOCSEARCH {

    /**
     * Coordinate descent for box constrained problems with unequal dynamically adjacent
     * steps along directions
     */
    template <typename FT> class VarCoorDesc : public COMPI::Solver<FT> {
    public:
        
        /**
         * Determines stopping conditions
         * @param xdiff - distance between next and previous x
         * @param fdiff - difference between next and previous f value
         * @param gran - current granularity vector
         * @param n - current step number
         */
        typedef std::function<bool(FT xdiff, FT fdiff, const std::vector<FT>& gran, FT fval, int n) > Stopper;

        /**
         * Watches the current step
         * @param xdiff - distance between next and previous x
         * @param fdiff - difference between next and previous f value
         * @param gran - current granularity vector
         * @param n - current step number
         */
        typedef std::function<void(FT xdiff, FT fdiff, const std::vector<FT>& gran, FT fval, int n) > Watcher;

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
         * @param stopper - reference to the stopper
         * @param ls - pointer to the line search
         */
        VarCoorDesc(const COMPI::MPProblem<FT>& prob, const Stopper& stopper) :
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
            const snowgoose::Box<double>& box = *(mProblem.mBox);
            std::vector<FT> sft;
            sft.assign(n, mOptions.mHInit);

            // One step
            auto step = [&] () {
                FT xd = 0;
                for (int i = 0; i < n; i++) {
                    FT h = sft[i];
                    FT dh = h *  mOptions.mShifts[i];
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
                        FT t = h * mOptions.mInc;                        
                        t = SGMIN(t, mOptions.mHUB);
                        sft[i] = t;
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
                        FT t = h * mOptions.mDec;                        
                        t = SGMAX(t, mOptions.mHLB);
                        sft[i] = t;
                        x[i] = tmp;
                    } else {
                        FT t = h * mOptions.mInc;                        
                        t = SGMIN(t, mOptions.mHUB);
                        sft[i] = t;
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
                } else {
                    FT H = snowgoose::VecUtils::maxAbs(n, sft.data(), nullptr);
                    if (H <= mOptions.mHLB)
                        break;                    
                }
                for(auto w : mWatchers) {
                    w(xdiff, fdiff, sft, fcur, sn);
                }
                if (mStopper(xdiff, fdiff, sft, fcur, sn)) {
                    break;
                }

            }
            v = fcur;
            return rv;
        }

        std::string about() const {
            std::ostringstream os;
            os << "Coordinate descent method with variable adaptation\n";
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

        Stopper mStopper;
        const COMPI::MPProblem<FT>& mProblem;
        LineSearch<FT>* mLS;
        Options mOptions;
        std::vector<Watcher> mWatchers;
        
    };
}

#endif /* GFSDESC_HPP */

