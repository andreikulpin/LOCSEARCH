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
         * @param stpn step number
         * @param gran - current granularity vector
         */
        using Watcher = std::function<void(FT fval, const FT* x, const std::vector<FT>& gran, int stpn) >;

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
             * Hooke-Jeeves multiplier (if <= 0 then don't do Hooke-Jeeves step)
             */
            FT mHJ = -1;
        };

        /**
         * The constructor
         * @param prob - reference to the problem
         */
        VarCoorDesc(const COMPI::MPProblem<FT>& prob) :
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
            const snowgoose::Box<double>& box = *(mProblem.mBox);
            std::vector<FT> sft;
            sft.assign(n, mOptions.mHInit);
            FT* xold = new FT[n];
            FT* xnew = new FT[n];

            auto inc = [this] (FT h) {
                FT t = h * mOptions.mInc;
                t = SGMIN(t, mOptions.mHUB);
                return t;
            };

            auto dec = [this](FT h) {
                FT t = h * mOptions.mDec;
                t = SGMAX(t, mOptions.mHLB);
                return t;
            };

            auto step = [&] () {
                int dir = 1;
                for (int i = 0; i < n;) {
                    FT dh = sft[i];
                    FT y = x[i] + dir * dh;
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
                            sft[i] = dec(dh);
                            i++;
                            dir = 1;
                        }
                    } else {
                        sft[i] = inc(dh);
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
                if (mOptions.mHJ > 0) {
                    snowgoose::VecUtils::vecCopy(n, x, xold);
                }
                step();
                if (fcur < fold) {
                    rv = true;
                    if (mOptions.mHJ > 0) {
                        FT l = mOptions.mHJ * snowgoose::VecUtils::maxAbs(n, sft.data(), nullptr);
                        while (true) {
                            for (int i = 0; i < n; i++) {
                                xnew[i] = x[i] + l * (x[i] - xold[i]);
                            }
                            snowgoose::BoxUtils::project(xnew, box);
                            FT fn = obj->func(xnew);
                            if (fn < fcur) {
                                snowgoose::VecUtils::vecCopy(n, xnew, x);
                                fcur = fn;
                                std::cout << "S: " << fn << ", l = " << l << "\n";
                                l *= 1.2;
                            } else {
                                std::cout << "F: " << fn << ", l = " << l << "\n";
                                break;
                                /*
                                if (l < mOptions.mHLB)
                                    break;
                                else
                                    l *= 0.5;
                                 */
                            }
                        }
                    }
                } else {
                    FT H = snowgoose::VecUtils::maxAbs(n, sft.data(), nullptr);
                    if (H <= mOptions.mHLB)
                        break;
                }
                for (auto w : mWatchers) {
                    w(fcur, x, sft, sn);
                }
                for (auto s : mStoppers) {
                    if (s(fcur, x, sn)) {
                        br = true;
                        break;
                    }
                }
            }
            v = fcur;
            delete [] xold;
            delete [] xnew;
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
        std::vector<Stopper> mStoppers;


    };
}

#endif /* GFSDESC_HPP */

