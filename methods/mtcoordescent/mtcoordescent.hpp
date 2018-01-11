/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   mtcoordescent.hpp
 * Author: mposypkin
 *
 * Created on October 26, 2017, 12:24 PM
 */

#ifndef MTCOORDESCENT_HPP
#define MTCOORDESCENT_HPP

#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <memory>
#include <omp.h>
#include <solver.hpp>
#include <common/dummyls.hpp>
#include <common/vec.hpp>
#include <box/boxutils.hpp>
#include <common/sgerrcheck.hpp>
#include <mpproblem.hpp>
#include <mputils.hpp>
#include <common/lineseach.hpp>


namespace LOCSEARCH {

    /**
     * Multithreaded advanced coordinate descent for box constrained problems with unequal dynamically adjacent
     * steps along directions combined 
     */
    template <typename FT> class MTCoordinateDescent : public COMPI::Solver<FT> {
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
         * Policy of vicinity adaptation
         */
        enum VicinityAdaptationPolicy {
            /**
             * No adaptation (vicinity is constant)
             */
            NO_ADAPTATION,
            /**
             * Adapting a whole vicinity is a box
             */
            UNIFORM_ADAPTATION,
            /**
             * Adaptation is along each axis
             */
            VARIABLE_ADAPTATION,
            /**
             * Assymetric adaptation (adaptation along each direction of search - two for each axis)
             */
            ASSYMETRIC_ADAPTATION
        };

        /**
         * Options for Gradient Box Descent method
         */
        struct Options {
            /**
             * Run in parallel mode or not
             */
            bool mParallelMode = true;
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
             * Lower bound on gradient
             */
            FT mGradLB = 1e-08;
            /**
             * Adaptation type
             */
            VicinityAdaptationPolicy mVicinityAdaptation = VicinityAdaptationPolicy::ASSYMETRIC_ADAPTATION;
        };

        /**
         * The constructor
         * @param prob - reference to the problem
         * @param stopper - reference to the stopper
         * @param ls - pointer to the line search
         */
        MTCoordinateDescent(const COMPI::MPProblem<FT>& prob) :
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
        bool search(FT* x, FT& v) override {
            bool rv = false;
            auto obj = mProblem.mObjectives.at(0);
            snowgoose::BoxUtils::project(x, *(mProblem.mBox));
            FT fcur = obj->func(x);
            const int n = mProblem.mVarTypes.size();
            const int n2 = 2 * n;
            const snowgoose::Box<double>& box = *(mProblem.mBox);
            std::vector<FT> sft(n2, mOptions.mHInit);
            std::vector<FT> xold(x, x + n);
            std::vector<FT> grad(n);
            std::vector< std::vector<FT> > xx(n2, xold);
            std::vector<FT> fv(n2);
            auto inc = [this] (FT & s) {
                s = std::min(s * mOptions.mInc, mOptions.mHUB);
            };
            auto dec = [this] (FT & s) {
                s = std::max(s * mOptions.mDec, mOptions.mHLB);
            };
            while (true) {

#pragma omp parallel for if(mOptions.mParallelMode)
                for (int i = 0; i < n2; i++) {
                    xx[i] = xold;
                    const int I = i / 2;
                    FT xi;
                    xi = xold[I] + (2 * (i % 2) - 1) * sft[i];
                    xi = std::min(xi, box.mB[I]);
                    xi = std::max(xi, box.mA[I]);
                    xx[i][I] = xi;
                    fv[i] = obj->func(xx[i].data());
                }
                // end parallel loop

                for (int i = 0; i < n; i++) {
                    const int j = 2 * i;
                    grad[i] = (fv[j + 1] - fv[j]) / (sft[j] + sft[j + 1]);
                }
                for (int i = 0; i < n2; i++) {
                    if (fv[i] < fcur) {
                        inc(sft[i]);
                    } else {
                        dec(sft[i]);
                    }
                }
                auto it = std::min_element(fv.begin(), fv.end());
                const FT fvn = *it;
                if (fvn < fcur) {
                    fcur = fvn;
                    const int I = std::distance(fv.begin(), it);
                    xold = xx[I];
                }
                auto itmax = std::max_element(sft.begin(), sft.end());
                if (*itmax <= mOptions.mHLB)
                    break;
                FT gnorm = 0;
                std::for_each(grad.begin(), grad.end(), [&gnorm](FT & el) {
                    gnorm += el * el;
                });
                if (gnorm * gnorm < mOptions.mGradLB)
                    break;
            }
            std::copy(xold.begin(), xold.end(), x);
            v = fcur;
            return rv;
        }

        std::string about() const {
            std::ostringstream os;
            os << "Advanced coordinate descent method with ";
            if (mOptions.mVicinityAdaptation == VicinityAdaptationPolicy::UNIFORM_ADAPTATION) {
                os << "uniform vicinity adaptation\n";
            } else if (mOptions.mVicinityAdaptation == VicinityAdaptationPolicy::VARIABLE_ADAPTATION) {
                os << "variable vicinity adaptation\n";
            } else {
                os << "no vicinity adaptation\n";
            }
            os << "Initial step = " << mOptions.mHInit << "\n";
            os << "Increment multiplier = " << mOptions.mInc << "\n";
            os << "Decrement multiplier = " << mOptions.mDec << "\n";
            os << "Upper bound on the vicinity size = " << mOptions.mHUB << "\n";
            os << "Lower bound on the vicinity size = " << mOptions.mHLB << "\n";
            os << "Lower bound on the gradient norm = " << mOptions.mGradLB << "\n";
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

        const COMPI::MPProblem<FT>& mProblem;
        Options mOptions;
        std::vector<Watcher> mWatchers;
        std::unique_ptr<LineSearch<FT>> mLS;

    };
}


#endif /* MTCOORDESCENT_HPP */

