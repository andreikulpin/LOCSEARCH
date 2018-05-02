/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   acd_bench.hpp
 * Author: andrei
 *
 * Created on May 2, 2018, 8:30 PM
 */

#ifndef ACD_BENCH_HPP
#define ACD_BENCH_HPP

#include <sstream>
#include <vector>
#include <functional>
#include <memory>
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
     * Advanced coordinate descent for box constrained problems with unequal dynamically adjacent
     * steps along directions combined with a descent along Hooke-Jeeves or anti-pseudo-gradient direction
     */
    template <typename FT> class AdvancedCoordinateDescentBench : public COMPI::Solver<FT> {
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
         * Search types
         */
        enum SearchTypes {
            /**
             * No search along descent direction if performed
             */
            NO_DESCENT,
            /**
             * Search along Hooke-Jeeves direction
             */
            HOOKE_JEEVES,
            /**
             * Search along anti-pseudo-grad
             */
            PSEUDO_GRAD
        };

        /**
         * Policy of vicinity adaptation
         */
        enum VicinityAdaptationPolicy {
            /**
             * No adaptation (vicinity is constant)
             */
            NO_ADAPTATION,
            /**
             * Vicinity is a box
             */
            UNIFORM_ADAPTATION,
            /**
             * Adaptation is along each direction
             */
            VARIABLE_ADAPTATION
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
            /**
             * Lower bound on gradient
             */
            FT mGradLB = 1e-08;
            /**
             * Searh type
             */
            SearchTypes mSearchType = SearchTypes::PSEUDO_GRAD;
            /**
             * Adaptation type
             */
            VicinityAdaptationPolicy mVicinityAdaptation = VicinityAdaptationPolicy::VARIABLE_ADAPTATION;
            /**
             * Trace on/off
             */
            bool mDoTracing = false;
            /**
             * Stop criterion
             */
            FT mEps = 0.6;
            /**
             * Max steps number
             */
            int maxStepsNumber = 100;
        };

        /**
         * The constructor
         * @param prob - reference to the problem
         * @param stopper - reference to the stopper
         * @param ls - pointer to the line search
         */
        AdvancedCoordinateDescentBench(const COMPI::MPProblem<FT>& prob, FT globMin) :
        mProblem(prob), mGlobMin(globMin) {
            unsigned int typ = COMPI::MPUtils::getProblemType(prob);
            SG_ASSERT(typ == COMPI::MPUtils::ProblemTypes::BOXCONSTR | COMPI::MPUtils::ProblemTypes::CONTINUOUS | COMPI::MPUtils::ProblemTypes::SINGLEOBJ);
        }

        /**
         * Retrieve the pointer to the line search
         * @return 
         */
        std::unique_ptr<LineSearch<FT>>&getLineSearch() {
            return mLS;
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
            int n = mProblem.mVarTypes.size();
            const snowgoose::Box<double>& box = *(mProblem.mBox);
            const int n2 = 2 * n;
            std::vector<FT> sft(n2, mOptions.mHInit);
            FT H = mOptions.mHInit;
            FT* xold = new FT[n];
            FT* grad = new FT[n];
            FT* ndir = new FT[n];

            auto inc = [this] (FT & s) {
                if (mOptions.mVicinityAdaptation == VicinityAdaptationPolicy::VARIABLE_ADAPTATION) {
                    s = std::min(s * mOptions.mInc, mOptions.mHUB);
                }
            };
            auto dec = [this] (FT & s) {
                if (mOptions.mVicinityAdaptation == VicinityAdaptationPolicy::VARIABLE_ADAPTATION) {
                    s = std::max(s * mOptions.mDec, mOptions.mHLB);
                }
            };

            auto step = [&] () {
                int dir = 1;
                for (int i = 0; i < n;) {
                    const int sfti = (dir > 0) ? 2 * i : 2 * i + 1;
                    const FT h = sft[sfti];
                    FT y = x[i] + dir * h;
                    y = std::max(y, box.mA[i]);
                    y = std::min(y, box.mB[i]);
                    const FT tmp = x[i];
                    x[i] = y;
                    const FT fn = obj->func(x);
                    const FT dx = y - tmp;
                    const FT ng = (dx == 0) ? 0 : (fn - fcur) / dx;
                    if (dir == 1) {
                        grad[i] = ng;
                    } else {
                        grad[i] = 0.5 * (grad[i] + ng);
                    }
                    if (fn >= fcur) {
                        x[i] = tmp;
                        dec(sft[sfti]);
                        if (dir == 1)
                            dir = -1;
                        else {
                            i++;
                            dir = 1;
                        }
                    } else {
                        inc(sft[sfti]);
                        fcur = fn;
                        dir = 1;
                        i++;
                    }
                }
            };

            int stepNum = 0;
            bool br = false;
            while (!br) {
                stepNum++;
                FT fold = fcur;
                if (mOptions.mSearchType == SearchTypes::HOOKE_JEEVES)
                    snowgoose::VecUtils::vecCopy(n, x, xold);
                step();
                if (mOptions.mVicinityAdaptation == VicinityAdaptationPolicy::UNIFORM_ADAPTATION) {
                    if (fcur < fold) {
                        H *= mOptions.mInc;
                        H = std::min(H, mOptions.mHUB);
                    } else {
                        H *= mOptions.mDec;
                        H = std::max(H, mOptions.mHLB);
                    }
                    sft.assign(n2, H);
                }
                FT gnorm = snowgoose::VecUtils::vecNormTwo(n, grad);
                if (mOptions.mDoTracing) {
                    std::cout << "After step fcur = " << fcur << ", fold = " << fold << "\n";
                    std::cout << "GNorm = " << gnorm << "\n";
                }
                if (gnorm <= mOptions.mGradLB) {
                    break;
                }
                if (mOptions.mDoTracing) {
                    std::cout << "Vicinity size = " << snowgoose::VecUtils::maxAbs(n2, sft.data(), nullptr) << "\n";
                }

                if (fcur == fold) {
                    FT vs = snowgoose::VecUtils::maxAbs(n2, sft.data(), nullptr);
                    if (vs <= mOptions.mHLB)
                        break;
                } else {
                    rv = true;
                }

                if (mOptions.mSearchType != SearchTypes::NO_DESCENT) {
                    if (mOptions.mSearchType == SearchTypes::PSEUDO_GRAD) {
                        SG_ASSERT(mLS);
                        snowgoose::VecUtils::vecMult(n, grad, -1., ndir);
                        if (mOptions.mDoTracing) {
                            std::cout << "Start Line search along anti pseudo grad\n";
                        }
                    } else if (mOptions.mSearchType == SearchTypes::HOOKE_JEEVES) {
                        if (fcur < fold) {
                            SG_ASSERT(mLS);
                            snowgoose::VecUtils::vecSaxpy(n, x, xold, -1., ndir);
                            if (mOptions.mDoTracing) {
                                std::cout << "Start Line search along Hooke-Jeeves direction\n";
                            }
                        }
                    }
                    snowgoose::BoxUtils::projectDirection(ndir, x, box);
                    if (mOptions.mDoTracing) {
                        std::cout << "x = " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
                        std::cout << "dir = " << snowgoose::VecUtils::vecPrint(n, ndir) << "\n";
                    }
                    mLS.get()->search(ndir, x, fcur);
                }

                if (stepNum >= mOptions.maxStepsNumber) {
                    br = true;
                    //std::cout << "Stopped as number of steps was too big\n";
                }
                
                if (SGABS(fcur - mGlobMin) < mOptions.mEps) {
                    br = true;
                    //std::cout << "Stopped as result reached target accuracy\n";
                }

                for (auto w : mWatchers) {
                    w(fcur, x, sft, stepNum);
                }
                for (auto s : mStoppers) {
                    if (s(fcur, x, stepNum)) {
                        br = true;
                        break;
                    }
                }
            }
            v = fcur;
            delete [] grad;
            delete [] xold;
            delete [] ndir;
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
            if (mOptions.mSearchType == SearchTypes::PSEUDO_GRAD) {
                SG_ASSERT(mLS);
                os << "Line search along anti pseudo-gradient direction is " << mLS.get()->about() << "\n";
            } else if (mOptions.mSearchType == SearchTypes::HOOKE_JEEVES) {
                SG_ASSERT(mLS);
                os << "Line search along Hooke-Jeeves direction is " << mLS.get()->about() << "\n";
            }
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
        std::vector<Stopper> mStoppers;
        std::vector<Watcher> mWatchers;
        std::unique_ptr<LineSearch < FT>> mLS;
        FT mGlobMin;
    };
}

#endif /* ACD_BENCH_HPP */

